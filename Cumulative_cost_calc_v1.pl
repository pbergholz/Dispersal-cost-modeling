#!/usr/bin/perl -w

#This script is to calculate cumulative movement costs amongst field sampling sites for bacteria transported by a variety of hypothetical vehicles.
#This script assumes that you have already run a cost map and attraction map calculation script like Cost_map_calculator.v3.pl.

#The script requires these inputs:
# INPUT 1: (-m) A vector map of field site points as an ASCII file.  
# INPUT 2: (-i) A text input file provided to Cost_map_calculator.pl.  Only the model names and a couple of the input maps are actually used, but the script assumes that the whole input file is there.
# OUTPUT 1: (-o) This is a file in which the work of the script will be summarized.
# OUTPUT COST TABLES:  Cumulative cost distance matrices will be reported as one file for each map.
# OUTPUT DISPERSAL TABLES: These association matrices will be calculated by subtracting destination attraction map values from source load map values and dividing by cost distance estimates.

###NOTE: The dispersal values increase with increasing dispersal rate, but cost-distance increases with decreasing dispersal rate.  A Mantel test can only estimate the significance of positive correlations among association matrices (in vegan, at least) so the measurement of actual dispersal must positively correlate with the model prediction (dispersal in the case of CPS).

###NOTE2:  This script only produces estimates from the multitude of landscape resistance models.  There is no functionality for the alternatives to landscape resistance (IBD and barrier).


#### MODULES ####

use Getopt::Std;
use v5.10;
use Data::Dumper;

#################

#### GET OPTIONS ####

getopt('imo');

if(!$opt_i)
	{
	die "No input file supplied.\n\nusage\:\t\tCumulative_cost_calc_v1.pl \-m \<Vector point map of field sites as ASCII file\> \-i \<LIST OF MAPS AND MODELS\> \-o \<OPTIONAL\: Model Generations Reports\>\n";
	}
if(!$opt_m)
	{
	die "No vector point map name supplied.\n\nusage\:\t\tCumulative_cost_calc_v1.pl \-m \<Vector point map of field sites\> \-i \<List of Maps and Models\> \-o \<OPTIONAL\: Model Generation Report\>\n";
	}
if(!$opt_o)
	{
	$opt_o="Cumcost_report\.out";  ### If output filename string is unsupplied, copy the input filename and tack '.out' to the end of it.
	}

#####################

#### PARSE INPUT FILE #####

open INPUT,"<$opt_i" or die "Can't open input file\!\!\n\n";

open OUTPUT, ">$opt_o" or die "Can't create output file\!\!\n\n";

my $starttime = localtime;
print OUTPUT "Cost model estimation script started $starttime\n";
print STDERR "Cost model estimation script started $starttime\n";

my $blockline='##############################';


my ($region) = <INPUT>=~/REGION\s+(\S+)/;
my ($sitesmap) = <INPUT>=~/SITES\s+(\S+)/;
my ($terrainmap) = <INPUT>=~/TERRAIN\s+(\S+)/;
my ($ibdmap) = <INPUT> =~/NULL\s+(\S+)/;
my ($sourcemap) = <INPUT>=~/SOURCES\s+(\S+)/;

#<INPUT>;  ###Throw out the COSTMAPS line as well.
#my (%baseline,%proximity,%barrier,%riparian,%density); ##Initialize global hashes associating map name with type.
#until ((my $record=<INPUT>)=~/ATTRACTMAPS/) ##Put the chunk in $record until it contains 'ATTRACTMAPS'
#	{
#	chomp $record;  ##Remove newline from end of chunk.
#	my ($maptype,$hashline)=split /\t/,$record; ###Break the lines down into type and commma-delimited hash initialization lines.
#	my @hash_line = split ',',$hashline; ###Split effect levels and mapnames into array units.
#		$hash_line=join ',',@hashline;
#	$maptype=lc $maptype;  ### Make $maptype all lower case.
#		print STDERR "$maptype\t$hash_line\n";
#	for (my $i=0;$i<scalar(@hash_line);$i+=2) ### Go through the array by pairs of effect and mapname
#		{
#		given ($maptype){   ###given is the new switch-like control structure built into Perl 5.10 and higher.
#		when('baseline'){$baseline{lc $hash_line[$i]}=$hash_line[$i+1];}
#		when('riparian'){$riparian{lc $hash_line[$i]}=$hash_line[$i+1];}
#		when('barrier'){$barrier{lc $hash_line[$i]}=$hash_line[$i+1];}
#		when('proximity'){$proximity{lc $hash_line[$i]}=$hash_line[$i+1];}
#		when('density'){$density{lc $hash_line[$i]}=$hash_line[$i+1];}
#		default {next;}} ### given requires 'use v5.10' in the module calls.
#		print OUTPUT "Logged input map: $hash_line[$i]\t$hash_line[$i+1]}\n";
#		}
#	}

#print STDERR "Extracted input map names...\n";
#print OUTPUT "Extracted input map names...\n";

until ((my $record=<INPUT>)=~/MODELS/)  ###Skip the section of the input file that lists attraction maps. ##Untested control structure to work around addition of attraction maps for the attraction calculator.
	{next;}

###########################

#### Test map hash initialization ####
#### Hashes should have effect names as keys and GRASS map names as values
#### Works: 4APR2012
#my $blockline = "##########################################";
#print OUTPUT "Baseline maps\n";
#foreach(keys %baseline)
 # {print OUTPUT "$_\:$baseline{$_}\n";}
#print OUTPUT "$blockline\n";
#print OUTPUT "Proximity maps\n";
#foreach(keys %proximity)
 # {print OUTPUT "$_\:$proximity{$_}\n";}
#print OUTPUT "$blockline\n";
#print OUTPUT "Riparian maps\n";
#foreach(keys %riparian)
 # {print OUTPUT "$_\:$riparian{$_}\n";}
#print OUTPUT "$blockline\n";
#print OUTPUT "Barrier maps\n";
#foreach(keys %barrier)
 # {print OUTPUT "$_\:$barrier{$_}\n";}
#print OUTPUT "$blockline\n";
#print OUTPUT "Density maps\n";
#foreach(keys %density)
 # {print OUTPUT "$_\:$density{$_}\n";}
#print OUTPUT "$blockline\n";
#print Dumper(%baseline);
######################################

#### Call R, import maps, calculate proximity functions #####
#### Test block  ###RSPerl works, but there's no functional mechanic for creating R objects in the Perl-initiated R environment => RReferences is widely regarded as non-functional, and this is the mechanism for exporting R objects to Perl correctly.
#### Next try: use perl to automate GRASS r.mapcalc directly.  4APR2012
##############################################################

###GRASS Automation###  ###On the upside, GRASS manages memory for raster map operations very well, whereas R does not.
my $modelheaders = <INPUT>;
chomp $modelheaders;
my @modelcolumns = split /\t/, lc $modelheaders;  #Grab the column titles.
#@modelcolumns = lc @modelcolumns;
my %modelhash;  ##This will store hashes for all models and their components for later use.
shift @modelcolumns; ##Shift out the one for the "Model" column.

###Specify Computational Region Explicitly ######
system "g.region res=3\.0 region=$region";  #####Fill in the region resolution here, set the boundaries to match a cost model map.  In an ideal world, the user has already done that, but it is easy to overlook.

### 13JAN2014

#Need to confirm that cumulative cost matrix calculator works.  Test this with only the first two models.
#Need to add code for load and attraction sampling
#Need to add code for dispersal matrix calculator.

#############################
while(my $record=<INPUT>)
	{
	chomp $record;
	my ($model,@effects) = split /\t/,$record;  ##Model number will be the modelhash key referencing a hash containing effectname keys and values will be effect levels.
	my $i=0;  ###We're gonna keep track of array elements with this counter.
	foreach(@modelcolumns) 
		{
		$modelhash{$model}->{"$_"}=lc $effects[$i]; ##Store a hash of model components associated with each model.
#		print STDERR "$_\t$effects[$i]\n";
		$i++;
		}
	}
print OUTPUT "Extracted model characteristics...\n";
print STDERR "Extracted model characteristics...\n";
#### TEST HASH FOR MODELS  ######
#print OUTPUT Dumper(%modelhash);
print OUTPUT "$blockline\n";
#print STDERR Dumper(%modelhash);
####Passed.  There might be some extra space at the end of the last entry in the column name array, so need to check for that. 5APR2012.
#################################



##Cost distance among sites!!!!  This will cost probably 12-24 hours to write.

#open MAP, "$opt_m" or die "Can't open ASCII vector point map\n\n";
#my @asciiheaders = split /\t/, <MAP>;	#I don't think there are headers in a standard ascii map output; ##Basic ASCII MAP is Xcoor\tYcoor\tcat.
#my %costmatrix; ### This will record the cumulative cost distance values for keys representing every site.  The values will be for that site as a start site and all others as stop sites.

foreach (sort {$a <=> $b} keys %modelhash)
	{
	my $modelname=$_;
	my %componenthash=%{$modelhash{$modelname}};
	my $modelmap = "Model".$modelname;
	my $modelsites = "Model".$modelname."\_sites";
	open MAP, "$opt_m" or die "Can't open ASCII vector point map\n\n";  #Think this needs to be here, or loop will only successfully run once.
	system "g.copy vect=$sitesmap\,$modelsites";  ##There is no overwrite flag for this operation.  This map may need to be deleted manually for each run of the script on a single mapset.
	while (my $record = <MAP>)
		{
		chomp $record;
		my ($xcoord,$ycoord,$cat) = split /\t/, $record;## Get each vector cat in turn.
		###Might be best to import the vector points map from ascii file for each model.
		system "v.extract --overwrite --quiet input=$modelsites list=$cat type=point output=Start$cat";  #Need the vector map file and a category number from the ascii vector map row
		wait;
		my $cumcostmap = "$modelmap"."cost"."start"."$cat";
		my $max_costdistance = 0; #Need to fill in NULL cost distance values with the vehicle-wise maximum
		given($componenthash{'vehicle'})
			{
			when('lt')
				{
				system "r.walk -k elevation=$terrainmap friction=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites max_cost=100000 walk_coeff=0.72,4.0,1.5,-2 slope_factor=-0.4663";
				wait;
				$max_costdistance=100000;
				}  
			when('mb')
				{
				system "r.cost -k input=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites max_cost=100000";
				wait;
				$max_costdistance=100000;
				}
			when('lb')
				{
				system "r.cost -k input=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites max_cost=100000";
				wait;
				$max_costdistance=100000;
				}
			when('st')
				{
				system "r.walk -k elevation=$terrainmap friction=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites max_cost=20000 walk_coeff=1.2,1,1,-1 slope_factor=-0.4663";  ##14JAN2014 - Changed from 2000 to 8000 to provide for at least some linkage in high cost maps, otherwise, what's the point?
				wait;
				$max_costdistance=20000;
				}
			}
		system "v.db.addcol map=$modelsites columns=\"Site$cat double precision\""; ##Add a column to the field site vector db to record cost distances from each point as a start point for this model.
		wait;
		system "v.what.rast vector=$modelsites raster=$cumcostmap column=Site$cat"; ##Sample the cost map at all point locations.
		wait;
		system "v.db.update map=$modelsites column=\"Site$cat\" value=$max_costdistance where=\"Site$cat IS NULL\"";
		wait;
		system "g.remove rast=$cumcostmap"; ###Remove the cost map;
		}

	##Get load and attraction values for each site.
	system "v.db.addcol map=$modelsites columns=\"Loadi double precision,Attractj double precision\"";
	wait;
	system "v.what.rast vector=$modelsites raster=$sourcemap column=Loadi";
	wait;
	$modelmap=~s/^Model/Attract/;  ##Change model map name to attraction map name.
	system "v.what.rast vector=$modelsites raster=$modelmap column=Attractj";
	wait;
	system "db.out.ogr input=$modelsites dsn=$modelsites.csv format=CSV"; #Export the cost distance data for this model as an ascii map.
	wait;
#	system "g.remove vect=$modelsites"; 
	print OUTPUT "$blockline\n";
	my $completiontime=localtime;
	print OUTPUT "Cumulative cost estimates $modelsites created at $completiontime\n";
	print STDERR "Cumulative cost estimates $modelsites created at $completiontime\n";
	close MAP;
	# Estimate Relative Dispersal Among Sites
	open CALCIN, "$modelsites.csv" or die "Can't open $modelsites.csv\n\n";
	open COSTMAT, ">$modelsites.cost.matrix.csv";
	open DISPMAT, ">$modelsites.dispersal.matrix.csv";
	my $headerstring = <CALCIN>;
	chomp $headerstring;
	my @headers = split /\,/, $headerstring;
	shift @headers;  #Discard the first two headers.
	shift @headers;
	my %calcinhash; 
	while (my $record = <CALCIN>)
		{
		chomp $record;
		my @siteline = split /\,/,$record;
		my $cat = shift @siteline;
		my $startsite = shift @siteline;
		$calcinhash{$cat} -> {'SITENAME'} = $startsite;
		for my $i (0..$#headers)
			{
			$calcinhash{$cat} -> {$headers[$i]} = $siteline[$i];  #Get values associated with Site cost-distance, Loadi and Attractj. 
			}
		}
#	print OUTPUT Dumper(%calcinhash);
	foreach (sort {$a cmp $b} keys %calcinhash)
		{
		my $startsite = $calcinhash{$_} -> {'SITENAME'};
		print COSTMAT "$startsite";
		print DISPMAT "$startsite";
		my %catline = %{$calcinhash{$_}};
		print OUTPUT Dumper(%catline);
		foreach (sort {$a cmp $b} keys %catline)
			{
			if ($_ eq 'Loadi' or $_ eq 'Attractj' or $_ eq 'SITENAME' or $_ eq 'Descriptio')
				{next;}
			my $loadi = $catline{'Loadi'};
			if ($loadi == 0)
				{$loadi+=0.1}
			my $cost_distance=$catline{$_}+0.1;
			my $stopsite = $_;
			$stopsite=~s/Site//;
			print COSTMAT "\,$catline{$_}";
			my $attractj = $calcinhash{$stopsite} -> {'Attractj'};
			if ($attractj==0)
				{$attractj+=0.1;}
			my $dispersal = ($loadi*$attractj)/$cost_distance;
			print DISPMAT "\,$dispersal";
			}
		print COSTMAT "\n";
		print DISPMAT "\n";
		}
	##################################
	}
my $completiontime = localtime;
print OUTPUT "Dispersal estimates script completed at $completiontime\n";
print STDERR "Dispersal estimates script completed at $completiontime\n";
###END OF SCRIPT ####

