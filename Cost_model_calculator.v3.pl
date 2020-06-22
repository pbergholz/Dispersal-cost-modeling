#!/usr/bin/perl -w

###This script is a module that calculates landscape resistance models.
###It should be considered one component of a connectivity calculator.
###This script only calculates landscape resistance models.
###For gravity models, habitat suitability, etc. other scripts will be needed.

## This script will organize cost model calculation for the CPS project.
## The script will send a series of commands for map calculation to a GRASS GIS shell.
## The key consideration is storage memory, therefore computational region for the GRASS location.
## Most regions I have considered so far are about 10 km on a side.
## So, 3mx3m resolution should bring diskspace for each map down to about 10-15 MB
## That's a workable size.
## This script will read a tab-delimited text file that dictates the essential
## model characteristics.  Those characteristics represent a set of cases for
## model calculation. The script will require that all baseline, addend and
## factor maps exist in GRASS. It will only run those through r.mapcalc to
## produce final cost maps.  It will feed the resulting maps to r.covar to
## produce a covariation matrix for all of the maps (I think there will be enough memory
## for that). The covariance matrix can then be used to understand what the major categories of
## models are from a quantitative standpoint.

## The calculation of cumulative cost distances by r.walk will be integrated into this script. !!! Checked with author of PATHMATRIX to see if there is a standalone version of that, there is, but it doesn't integrate DEM data, so I'm back to r.walk.

##My r.walk scripting block requires an ascii vector point map containing points representing sampled field centroids.  The script will extract each point, in turn, calculate r.walk from the extracted point as a starting point to all possible end points.  It will then write the data out to the ascii map file as a new column with a header indicating the start point and the values as extracts cumulative cost distance to each end point.


##The input file for this script is a tab-delimited file with three sections
## Section one should be titles in the first column: COSTMAPS (in all capitals)
## It should have one line for each needed map type:
##   Baseline, Riparian, Proximity, and Barrier
## The lines should contain:
##   type name<tab><comma-delimited list of interspersed effects and map names in type>
##   This format is like a comma-delimited hash initialization.

##Section two should be titled in the first column: ATTRACTMAPS (in all capitals)
## It should have one line for each needed map type:
##  Baseline, Area, Proximity, Habcover, Badcover
## The lines should be formatted the same as in section one:
##	map type name<tab><comma-delimited list of interspersed effects and map names in type>
##	The format is used as a comma-delimited hash initialization line.

## Section three should be titled in the first column: MODELS (in all capitals) 
## Section three is a tab delimited table one model description per line.
## The input file column titles are: Model (int),Vehicle (fact),
##   Riparian(fact),Proximity(fact),Barrier(fact),Attraction (fact).
## Attraction is not actually used in this calculator, but in the attraction map calculator. 
##   Levels for various factors are:
#BASELINE
## LT - Long dispersing terrestrial animals
## ST - Locally dispersing terrestrial animals
## LB - Birds that disperse near the ground, where landcover will impact dispersal by cost distance.
## MB - Birds that disperse long distances and are primarily limited by ground cover attraction rather than cost.
#RIPARIAN
## Adjacent - All lands adjacent to riparian zone cut dispersal in half regardless of their extent.
## DistBuffer - Lands adjacent to riparian out to 100 m cut dispersal cost in half.
## Topbuffer - Lands adjacent to riparian along the 50 yr flood height channel are half cost.
## None - No riparian effect.
#PROXIMITY
## High - Strong proximity effects are calculated by the script and used for the cost model.
## Low - Weak proximity effects are calculated by the script and used for the cost model.
## None - No proximity effects included.
#BARRIER
## Por - Porous barriers are included with each cell of barrier reset to cost = 200.
## Abs - Absolute barriers are included with each cell of barrier reset to cost = the max possible cost for the organism.
## None - Barriers have no effect.
#ATTRACTION
## Hab - % land cover from wildlife habitat has twice the attraction coefficient
## Area - Patch area has no effect.
## Norm - No changes to effects.
## None - No attraction effect, only cost distance function and load.

##The map files needed to run calculations are: Landuse baselines based on 
##   dispersal vehicle,
##Three different riparian effect maps.
##Proximity maps for:
## urban development and major roads, water, forests, pastures, wetlands.
##And Barrier maps for waterways and heavy traffic


#### MODULES ####

use Getopt::Std;
use v5.10;
use Data::Dumper;

#################

#### GET OPTIONS ####

getopt('mio');

if(!$opt_i)
	{
	die "No input file supplied.\n\nusage\:\t\tCost_model_calculator.v2.pl \-m \<Vector point map of field sites as ASCII file\> \-i \<LIST OF MAPS AND MODELS\> \-o \<OPTIONAL\: Model Generations Reports\>\n";
	}
if(!$opt_m)
	{
	die "No vector point map name supplied.\n\nusage\:\t\tCost_model_calculator.v2.pl \-m \<Vector point map of field sites\> \-i \<List of Maps and Models\> \-o \<OPTIONAL\: Model Generation Report\>\n";
	}
if(!$opt_o)
	{
	$opt_o=$opt_i."\.out";  ### If output filename string is unsupplied, copy the input filename and tack '.out' to the end of it.
	}

#####################

#### PARSE INPUT FILE #####

open INPUT,"<$opt_i" or die "Can't open input file\!\!\n\n";

open OUTPUT, ">$opt_o" or die "Can't create output file\!\!\n\n";

my $starttime = localtime;
print OUTPUT "Cost model estimation script started $starttime\n";
print STDERR "Cost model estimation script started $starttime\n";

my ($sitesmap) = <INPUT>=~/SITES\s+(\S+)/;
my ($terrainmap) = <INPUT>=~/TERRAIN\s+(\S+)/;

<INPUT>;  ###We don't need the null cost (distance) map for this script.
<INPUT>;  ###Throw out the COSTMAPS line as well.
my (%baseline,%proximity,%barrier,%riparian,%density); ##Initialize global hashes associating map name with type.
until ((my $record=<INPUT>)=~/ATTRACTMAPS/) ##Put the chunk in $record until it contains 'ATTRACTMAPS'
	{
	chomp $record;  ##Remove newline from end of chunk.
	my ($maptype,$hashline)=split /\t/,$record; ###Break the lines down into type and commma-delimited hash initialization lines.
	my @hash_line = split ',',$hashline; ###Split effect levels and mapnames into array units.
#		$hash_line=join ',',@hashline;
	$maptype=lc $maptype;  ### Make $maptype all lower case.
#		print STDERR "$maptype\t$hash_line\n";
	for (my $i=0;$i<scalar(@hash_line);$i+=2) ### Go through the array by pairs of effect and mapname
		{
		given ($maptype){   ###given is the new switch-like control structure built into Perl 5.10 and higher.
		when('baseline'){$baseline{lc $hash_line[$i]}=$hash_line[$i+1];}
		when('riparian'){$riparian{lc $hash_line[$i]}=$hash_line[$i+1];}
		when('barrier'){$barrier{lc $hash_line[$i]}=$hash_line[$i+1];}
		when('proximity'){$proximity{lc $hash_line[$i]}=$hash_line[$i+1];}
		when('density'){$density{lc $hash_line[$i]}=$hash_line[$i+1];}
		default {next;}} ### given requires 'use v5.10' in the module calls.
#		print OUTPUT "Logged input map: $hash_line[$i]\t$hash_line[$i+1]}\n";
		}
	}

print STDERR "Extracted input map names...\n";
print OUTPUT "Extracted input map names...\n";

until ((my $record=<INPUT>)=~/MODELS/)  ###Skip the section of the input file that lists attraction maps. ##Untested control structure to work around addition of attraction maps for the attraction calculator.
	{next;}

###########################

#### Test map hash initialization ####
#### Hashes should have effect names as keys and GRASS map names as values
#### Works: 4APR2012
my $blockline = "##########################################";
print OUTPUT "Baseline maps\n";
foreach(keys %baseline)
  {print OUTPUT "$_\:$baseline{$_}\n";}
print OUTPUT "$blockline\n";
print OUTPUT "Proximity maps\n";
foreach(keys %proximity)
  {print OUTPUT "$_\:$proximity{$_}\n";}
print OUTPUT "$blockline\n";
print OUTPUT "Riparian maps\n";
foreach(keys %riparian)
  {print OUTPUT "$_\:$riparian{$_}\n";}
print OUTPUT "$blockline\n";
print OUTPUT "Barrier maps\n";
foreach(keys %barrier)
  {print OUTPUT "$_\:$barrier{$_}\n";}
print OUTPUT "$blockline\n";
print OUTPUT "Density maps\n";
foreach(keys %density)
  {print OUTPUT "$_\:$density{$_}\n";}
print OUTPUT "$blockline\n";
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
system "g.region res=3\.0 rast=$baseline{'lt'}";  #####Fill in the region resolution here, set the boundaries to match one of the input raster maps.  In an ideal world, the user has already done that, but it is easy to over look.


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
print OUTPUT Dumper(%modelhash);
print OUTPUT "$blockline\n";
#print STDERR Dumper(%modelhash);
####Passed.  There might be some extra space at the end of the last entry in the column name array, so need to check for that. 5APR2012.
#################################


###Proximity effect calculation###

foreach(keys %proximity)
	{

	given(lc $_)
		{
		when('urban')
			{
			print STDERR "Calculating urban edge effect...\n";
			print OUTPUT "Calculating urban edge effect...\n";
			system "r.mapcalc \"LT_low_urban_edge_effect=0-(-10*(exp(0.01*250)/(exp(0.01*$proximity{$_})+exp(0.01*250)-1)))\"";###Low...confirmed function shape 31JUL2012
#			wait;
			system "r.mapcalc \"LT_high_urban_edge_effect=0-(-10*(exp(0.075*500)/(exp(0.075*$proximity{$_})+exp(0.075*500)-1)))\"";###High...confirmed function shape 31JUL2012			
			wait;
			}
		when('woody') ##I had noted that the effects of proximity to scrubland and woody wetland should be similar for LT animals.  I think that rather than calculating the same function on each land cover type in succession, I should generate a raster map that contains all three.
			{
			print STDERR "Calculating forest edge effect...\n";
			system "r.mapcalc \"LT_low_forest_edge_effect=10*(-1+exp(0.02*350)+exp(0.02*$proximity{$_}))/((exp(0.02*$proximity{$_})+1*exp(0.02*350)*10))\"";#Low...confirmed function shape 31JUL2012
#			wait;
			system "r.mapcalc \"LT_high_forest_edge_effect=10*((exp(0.02*$proximity{$_}))/(1+exp(0.02*$proximity{$_})))\"";#High...confirmed function shape 31JUL2012
			wait;
			}
		when('water')
			{
			print STDERR "Calculating water edge effect...\n";
			system "r.mapcalc \"LT_low_water_edge_effect=5*(-1+exp(0.005*250)+exp(0.005*$proximity{$_}))/(exp(0.005*$proximity{$_})-1+exp(0.005*250)*5)\""; ##Low...confirmed function shape 31JUL2012
#			wait;
			system "r.mapcalc \"LT_high_water_edge_effect=5*((exp(0.005*$proximity{$_}))/(1+exp(0.005*$proximity{$_})))\"";##High...confirmed function shape 31JUL2012
			wait;
			}
		when('rodentcover')  ###Rodent cover includes pastures, grasslands, scrublands, classes 1-4 of roads
			{
			print STDERR "Calculating habitat edge effects for small mammals...\n";
			system "r.mapcalc \"ST_low_cover_edge_effect=20*(-1+exp(0.075*50)+exp(0.075*$proximity{$_}))/((exp(0.075*$proximity{$_})+1*exp(0.075*50)*20))\""; ###Low...confirmed function shape 31JUL2012
#			wait;
			system "r.mapcalc \"ST_high_cover_edge_effect=20*((exp(0.025*$proximity{$_}))/(1+exp(0.025*$proximity{$_})))\""; ###High...confirmed function shape 31JUL2012
			wait;
			}
		when('lbroost')
			{print STDERR "Calculating proximity edge effects for flocking avian granivore/insectivores...\n";
			system "r.mapcalc \"LB_low_cover_edge_effect=2.0*(-1.0+exp(0.008*500)+exp(0.008*$proximity{$_}))/(exp(0.008*$proximity{$_})-10+exp(0.008*500)*5)\"";  ##Low...confirmed function shape 13MAR2013
#			wait;
			system "r.mapcalc \"LB_high_cover_edge_effect=5.0*((exp(0.01*$proximity{$_}))/(10+exp(0.01*$proximity{$_})))\""; ##High...confirmed function shape 13MAR2013.
			wait;
			}
		default{print STDERR "$_\:\nProximity function does not exist\n\n";print OUTPUT "ERROR\;: $_\:\nProximity function does not exist\n";next;}
		}
	}

print STDERR "Calculated all proximity effects successfully!\n";
print OUTPUT "Calculated all proximity effects successfully!\n";
print OUTPUT "$blockline\n";
###Percent land cover effects calculation######
####LT only respond to landcover in attraction models but all others factor it into cost.  This is supposed to reflect different scales of landscape observation.  ST have more pathways available in high habitat coverage, LB and MB prefer to move through these areas because they can stop over during foraging flight to rest and eat on their way to better locations.

foreach (keys %density)
	{
	given (lc $_)
		{
		print STDERR "Calculating density effects for movement costs...\n";
		print OUTPUT "Calculating density effects for movement costs...\n";
		when('mb')
			{
			##Rescore nulls to high cost!
#			print STDERR "Calculating habitat density effect on cost...\n";
			system "r.mapcalc \"MB_low_hab_density_effect=0-(-3*(exp(0.1*35)/(exp(0.1*$density{$_})+(exp(0.1*35)-3))))\""; ###Low...confirmed function shape 31JUL2012
			wait;
			system "r.mapcalc \"MB_low_hab_density_effect=if(isnull(MB_low_hab_density_effect),5,MB_low_hab_density_effect)\"";
			wait;
			system "r.mapcalc \"MB_high_hab_density_effect=0-(-3*exp(0.1*80)/(exp(0.1*$density{$_})+exp(0.1*80)-3))\"";  ###High...confirmed function shape 31JUL2012
			wait;
			system "r.mapcalc \"MB_high_hab_density_effect=if(isnull(MB_high_hab_density_effect),10,MB_high_hab_density_effect)\"";
			wait;
			}
		when('st')
			{
			system "r.mapcalc \"ST_low_hab_density_effect=0-(-10*exp(0.1*25)/(exp(0.1*$density{$_})+(exp(0.1*25))))\"";  ##Low...confirmed function shape 3AUG2012
			wait;
			system "r.mapcalc \"ST_low_hab_density_effect=if(isnull(ST_low_hab_density_effect),10,ST_low_hab_density_effect)\"";
			wait;
			system "r.mapcalc \"ST_high_hab_density_effect=0-(-10*exp(0.30*70)/(exp(0.30*$density{$_})+(exp(0.30*70))))\"";  ##High...confirmed function shape 3AUG2012
			wait;
			system "r.mapcalc \"ST_high_hab_density_effect=if(isnull(ST_high_hab_density_effect),20,ST_high_hab_density_effect)\"";
			wait;
			}
		when('lb')
			{
			system "r.mapcalc \"LB_low_density_effect=0-(-5.0*(exp(0.10*5))/(exp(0.10*($density{$_}))+(exp(0.10*5)-1.0)))\"";  ##Low...confirmed function shape 13MAR2013.
			wait;
			system "r.mapcalc \"LB_low_density_effect=if(isnull(LB_low_density_effect),5,LB_low_density_effect)\"";
			wait;
			system "r.mapcalc \"LB_high_density_effect=0-(-5.0*(exp(0.15*50))/(exp(0.15*($density{$_}))+(exp(0.15*50)-1.0)))\""; ##High..confirmed function shape 13MAR2013.
			wait;
			system "r.mapcalc \"LB_high_density_effect=if(isnull(LB_low_density_effect),10,LB_high_density_effect)\"";
			wait;
			}
		default{print STDERR "$_\:\nDensity function does not exist\n\n";print OUTPUT "ERROR\: $_\:\nDensity function does not exist\n";next;}	
		}
	}
print OUTPUT "Successfully calculated density effects on movement costs...\n";
print OUTPUT "$blockline\n";
##################################

### Substitute 0 for NULL in all partial maps ####
### Conditionals result in NULL if x is NULL, so doing math with conditions in a partial map will erase the basemap.

### This block assumes that only the barrier and the riparian maps have NULL values.
print OUTPUT "Recoding riparian maps...\n";

foreach (keys %riparian)
	{
	my $mapout=$riparian{$_};
#	print STDERR "Input name\:\t$riparian{$_}\n";
#	print STDERR "Held\:\t$mapout\n";
	$mapout = "$mapout"."recode";
	system "r.mapcalc \"$mapout=if\(isnull\($riparian{$_}\),0,$riparian{$_}\)\""; 
	wait;
	$riparian{$_}=$mapout; ##Replaces original riparian map with new riparian map, only in cases where original map is from a different mapset.
	}
print OUTPUT "All riparian maps recoded for 0 in place of NULL\n";
print OUTPUT "Recoding barrier maps...\n";
foreach (keys %barrier)
	{
	my $mapout = $barrier{$_};
#	print STDERR "Input name\:\t$barrier{$_}\n";
#	print STDERR "Held\:\t$mapout\n";
	$mapout="$mapout"."recode";
	system "r.mapcalc \"$mapout=if\(isnull\($barrier{$_}\),0,$barrier{$_}\)\""; 
	wait;
	$barrier{$_}=$mapout; ##See comment line 206.
	print STDERR "Barrier map $mapout recoded...\n";
	print OUTPUT "Barrier map $mapout recoded...\n";
	}
print OUTPUT "All barrier maps recoded for 0 in place of NULL\n";
print OUTPUT "$blockline\n";
##################################################

### Cost Model Calculation ####  ###The order of operations here matters a lot.
### This faithfully reproduces the cost map calculated in R for the LT (r=0.99994). Subtraction of the model calculated here, from the one that I calculated by hand produce differences ranging from 2x10^-6 to 1x10^-13.

print OUTPUT "Beginning cost model calculations...\n";

foreach(sort {$a <=> $b} keys %modelhash)
	{
#	my $mapcalc_call;
	my $modelname = $_;
#	my $modelbar=$modelname."B";
#	my $modelrip=$modelname."R";
#	my $modelprox=$modelname."P";
	my %componenthash=%{$modelhash{$modelname}}; #Pop out the components for this model.
	print OUTPUT "$blockline\n";
	print OUTPUT "Calculating model$modelname...\n";
	given($componenthash{'vehicle'}){
		when('lt') ##(l)ong distance dispersing, (t)errestrial mammal similar to deer or feral swine.
			{
			print OUTPUT "Model $modelname\:\tCalculating base cost -- Long dispersing terrestrial...\n";
			system "r.mapcalc \"Model$modelname=$baseline{'lt'}\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=2,20,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=3,10,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=4,5,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=5,3,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'water'}\=\=2,3,Model$modelname\)\"";  ### defines baseline
			wait;
##need to figure out conditionals in r.mapcalc before proceeding. 
###if(x,a) is NULL if x is NULL, a if x is non-zero,0 otherwise
###if(x,a,b) is NULL if x is NULL, a if x is non-zero, b otherwise
###This is important.  I added a block of code above that converts the feature-based rasters from NULL to 0.
			given($componenthash{'proximity'}){
				when('high')
					{
					print OUTPUT "Model$modelname\:\tAdding proximity effects -- High cost...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+LT_high_forest_edge_effect+LT_high_urban_edge_effect+LT_high_water_edge_effect\"";
					wait;
					}
				when('low')
					{
					print OUTPUT "Model$modelname\:\tAdding proximity effects -- Low cost...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+LT_low_forest_edge_effect+LT_low_urban_edge_effect+LT_low_water_edge_effect\"";
					wait;
					}
				when('none')
					{next;}
				default{print OUTPUT "Proximity cost level doesn't exist\:\tModel$modelname\n";next;}
				}
			given($componenthash{'riparian'}){
				when('adjacent')
					{
					print OUTPUT "Model $modelname\:\tMultiplying in riparian effect -- Adjacent...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'adjacent'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('distbuffer')
					{
					print OUTPUT "Model$modelname\:\tMultiplying in riparian effect -- Distance-based buffer...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'distbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('topbuffer')
					{
					print OUTPUT "Model$modelname\:\tMultiplying in riparian effect -- Topo-based buffer...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'topbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('none')
					{next;}
				default{print OUTPUT "Riparian cost level doesn't exist\:\tModel$modelname\n";next;}			
					}
			given($componenthash{'barrier'}){
				when('por')
					{
					print OUTPUT "Model $modelname\:\tAdding in barriers -- Porous...\n";
					system "r.mapcalc \"Model$modelname=if\($barrier{'water'}\=\=1,200,Model$modelname\)\"";
					wait;
					system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=1,200,Model$modelname\)\"";			
					wait;
					}
				when('abs')
					{
					print OUTPUT "Model $modelname\:\tAdding in barriers -- Absolute...\n";
					system "r.mapcalc \"Model$modelname=if\($barrier{'water'}\=\=1,40000,Model$modelname\)\"";
					wait;
					system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=1,40000,Model$modelname\)\"";
					wait;
					}
				when('none')
					{next;}
				default{print OUTPUT "Barrier level doesn't exist\:\tModel$modelname\n"; next;} ##Probably ought to have these printed to a file specific for critical error messages.
				}
			}
		when('mb')
			{
			print OUTPUT "Model $modelname\:\tCalculating base cost -- Migratory birds...\n\n";
			system "r.mapcalc \"Model$modelname=$baseline{'mb'}\"";
			given($componenthash{'proximity'})
				{
				when('high')
					{
					print OUTPUT "Adding in land cover density effects -- High...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+MB_high_hab_density_effect\"";
					wait;
					}
				when('low')
					{
					print OUTPUT "Adding in land cover density effects -- Low...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+MB_low_hab_density_effect\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Cost land cover density effect level doesn't exist\:\tModel$modelname\n";next;};
				}
			given($componenthash{'riparian'})
				{
				when('adjacent')
					{
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'adjacent'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('distbuffer')
					{
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'distbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('topbuffer')
					{
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'topbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Cost riparian effect level doesn't exist\:\tModel$modelname\n"; next;}
				}
			}
############################  STILL NEEDS WRITING: 15JAN2013 ################
		when('lb')  
			{
			print OUTPUT "Calculating baseline cost map -- Flocking granivores and insectivores\n\n";
			system "r.mapcalc \"Model$modelname=$baseline{'lb'}\"";
			wait;
			given($componenthash{'proximity'}){
				when('high')
					{
					print OUTPUT "Model$modelname\:\tAdding proximity effects -- High cost...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+LB_high_cover_edge_effect+LB_high_density_effect\"";
					wait;
					}
				when('low')
					{
					print OUTPUT "Model$modelname\:\tAdding proximity effects -- Low cost...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+LB_low_cover_edge_effect+LB_low_density_effect\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Proximity cost level doesn't exist\:\tModel$modelname\n";next;}
				}
			given($componenthash{'riparian'}){
				when('adjacent')
					{
					print OUTPUT "Model $modelname\:\tMultiplying in riparian effect -- Adjacent...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'adjacent'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('distbuffer')
					{
					print OUTPUT "Model$modelname\:\tMultiplying in riparian effect -- Distance-based buffer...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'distbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('topbuffer')
					{
					print OUTPUT "Model$modelname\:\tMultiplying in riparian effect -- Topo-based buffer...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'topbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Riparian cost level doesn't exist\:\tModel$modelname\n";next;}			
				}						
			}
#############################################################################
		when('st')
			{
			print OUTPUT "Calculating baseline cost map -- Short Dispersing Terrestrial\n\n";
			system "r.mapcalc \"Model$modelname=$baseline{'st'}\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=2,20,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=3,10,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=4,10,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=5,10,Model$modelname\)\"";
			wait;
			system "r.mapcalc \"Model$modelname=if\($barrier{'water'}\=\=2,40,Model$modelname\)\"";  ### defines baseline
			wait;
			given($componenthash{'proximity'}){
				when('high')
					{
					print OUTPUT "Model$modelname\:\tAdding proximity effects -- High cost...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+ST_high_cover_edge_effect+ST_high_hab_density_effect\"";
					wait;
					}
				when('low')
					{
					print OUTPUT "Model$modelname\:\tAdding proximity effects -- Low cost...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname+ST_low_cover_edge_effect+ST_low_hab_density_effect\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Proximity cost level doesn't exist\:\tModel$modelname\n";next;}
				}
			given($componenthash{'riparian'}){
				when('adjacent')
					{
					print OUTPUT "Model $modelname\:\tMultiplying in riparian effect -- Adjacent...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'adjacent'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('distbuffer')
					{
					print OUTPUT "Model$modelname\:\tMultiplying in riparian effect -- Distance-based buffer...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'distbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('topbuffer')
					{
					print OUTPUT "Model$modelname\:\tMultiplying in riparian effect -- Topo-based buffer...\n";
					system "r.mapcalc \"Model$modelname=Model$modelname*if\($riparian{'topbuffer'}\=\=1,0.5,1\)\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Riparian cost level doesn't exist\:\tModel$modelname\n";next;}			
					}
			given($componenthash{'barrier'}){
				when('por')
					{
					print OUTPUT "Model $modelname\:\tAdding in barriers -- Porous...\n";
					system "r.mapcalc \"Model$modelname=if\($barrier{'water'}\=\=1,50,Model$modelname\)\"";
					wait;
					system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=1,50,Model$modelname\)\"";			
					wait;
					}
				when('abs')
					{
					print OUTPUT "Model $modelname\:\tAdding in barriers -- Absolute...\n";
					system "r.mapcalc \"Model$modelname=if\($barrier{'water'}\=\=1,1000,Model$modelname\)\"";
					wait;
					system "r.mapcalc \"Model$modelname=if\($barrier{'roads'}\=\=1,1000,Model$modelname\)\"";
					wait;
					}
				when('none')
					{next;}
				default{print STDERR "Barrier level doesn't exist\:\tModel$modelname\n"; next;} ##Probably ought to have these printed to a file specific for critical error messages.
				}
			}
		}
	}  ####What's going on here?  SOLVED 19MAR2013


##### Split the below out to a separate program named Cumulative_cost_calc.pl 22MAR2013


###Since I've produced a functional calculator for cost surfaces, I'm going to forego the next block of code about using r.walk from each source location, and I'm going to move straight to attraction and load calculators.  I'll come back to cost distance estimation when we're closer on the other functions that go into flux.

##Cost distance among sites!!!!  This will cost probably 12-24 hours to write.

#open MAP, "$opt_m" or die "Can't open ASCII vector point map\n\n";
#my @asciiheaders = split /\t/, <MAP>;	#I don't think there are headers in a standard ascii map output; ##Basic ASCII MAP is Xcoor\tYcoor\tcat.
#my %costmatrix; ### This will record the cumulative cost distance values for keys representing every site.  The values will be for that site as a start site and all others as stop sites.
#foreach (sort {$a <=> $b} keys %modelhash)
#	{
#	my $modelname=$_;
#	my %componenthash=%{$modelhash{$modelname}};
#	my $modelmap = "Model".$modelname;
#	my $modelsites = "Model".$modelname."\_sites";
#	system "g.copy vect=$sitesmap\,$modelsites";  ##There is no overwrite flag for this operation.  This map may need to be deleted manually for each run of the script on a single mapset.
#	while (my $record = <MAP>)
#		{
#		chomp $record;
#		my ($xcoord,$ycoord,$cat) = split /\t/, $record;## Get each vector cat in turn.
		###Might be best to import the vector points map from ascii file for each model.
#		system "v.extract --overwrite input=$modelsites list=$cat type=point output=Start$cat";  #Need the vector map file and a category number from the ascii vector map row
#		wait;
#		my $cumcostmap = "$modelmap"."cost"."start"."$cat";
#		given($componenthash{'vehicle'})
#			{
#			when('lt')
#				{
#				system "r.walk -k elevation=$terrainmap friction=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites max_cost=40000 walk_coeff=0.72,4.0,1.5,-2 slope_factor=-0.4663";
#				wait;
#				}  ###r.walk parameters need to include slope costs.
#			when('mb')
#				{
#				system "r.walk -k elevation=$terrainmap friction=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites max_cost=40000 walk_coeff=0.1,1,1,10 slope_factor=3.73";
#				wait;
#				}
#			when('lb')
#				{
#				system "r.walk -k elevation=$terrainmap friction=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites walk_coeff=0.1,1,1,10 slope_factor=3.73";
#				wait;
#				}
#			when('st')
#				{
#				system "r.walk -k elevation=$terrainmap friction=$modelmap output=$cumcostmap start_points=Start$cat stop_points=$modelsites walk_coeff=1.2,1,1,-1 slope_factor=-0.4663";
#				wait;
#				}
#			}
#		system "v.db.addcol map=$modelsites columns=Site$cat double precision"; ##Add a column to the field site vector db to record cost distances from each point as a start point for this model.
#		wait;
#		system "v.what.rast vector=$modelsites raster=$cumcostmap column=Site$cat"; ##Sample the cost map at all point locations.
#		wait;
#		system "g.remove rast=$cumcostmap" ###Remove the cost map;
#		}
#	system "v.out.ascii input=$modelsites output=\.\/$modelsites\.txt"; ##Export the cost distance data for this model as an ascii map.
#	wait;
#	system "g.remove vect=$modelsites"; 
#	print "$blockline\n";
#	my $completiontime=localtime;
#	print OUTPUT "Cost estimates $modelsites created at $completiontime\n";
#	print STDERR "Cost estimates $modelsites created at $completiontime\n";
#	}
my $completiontime = localtime;
print OUTPUT "Cost estimates script completed at $completiontime\n";
print STDERR "Cost estimates script completed at $completiontime\n";
###END OF SCRIPT ####


