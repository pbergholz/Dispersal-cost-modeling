#!/usr/bin/perl -w

#This script is intended to function much like the cost models function.  It will take a file containing both a list of map names and a list dispersal vehicles (e.g., long dispersing terrestrial 'ruminant' class animals, short dispersing terrestrial 'rodent' class animals, near ground moving 'starling' avian class animals and high altitude dispersing 'goose' class avians).

##The script will calculate gravity functions that are weighted based on land cover area and density of favorable landcover near the target cropland.

### The inputs are largely based on papers about landscape factors that maximize probability of observation for species.  These are then parameterized based on data on the dispersal distances, nesting habitat and food preferences.  Croplands have a baseline attractiveness modified by patch area, but get no other direct modifications only indirect modifications based on their surroundings.

### INPUT MAPS ###

#Baseline attractiveness coefficient from landcover (species specific), scaled as 0.1 to 1.
#Patch area attractiveness coefficient (species specific function), scaled as 0.1 to 1.
#Habitat % land cover attractiveness coefficient (species specific function and habitat), scaled as 0.1 to 1.
#Repulsion % landcover attractiveness coefficient (species specific function and repulsive patches), scaled as 0.1 to 1.
#Exponential decay of attractiveness coefficient away from croplands (the focus of our study), scaled as 0.5 to 1.

##The maps have a minimum of 0.1, so that even in very unattractive areas, there is still some potential.  The middle of cropland fields would generally be scored as near zero for deer, though under my models, the edges can be very attractive. I want a little potential attractiveness within them.

##Call Modules

use v5.10;
use Getopt::Std;
use Data::Dumper;

### GET OPTIONS

getopt('io');

if(!$opt_i)
	{
	die "No input file supplied.\n\nusage\:\t\tattraction_calculator.pl \-i \<LIST OF MAPS AND MODELS\> \-o \<OPTIONAL\: Model Generations Reports\>\n";
	}

if(!$opt_o)  ###Output files should include at least an correlation coef matrix for all pairwise comparisons of all attraction models generated in the script. 
	{    ### The idea is to generate a PCA or PCoA of all of the maps to understand which classes of models are really different.
	$opt_o=$opt_i."\.out";  ### If output filename string is unsupplied, copy the input filename and tack '.out' to the end of it.
	}

#### PARSE INPUT FILE #####   ### This block was copied from the cost model calculator, so needs to be modified for this script.  12JUL2012.

open INPUT,"<$opt_i" or die "Can't open input file\!\!\n\n";

open OUTPUT, ">$opt_o" or die "Can't create output file\!\!\n\n";

my $starttime = localtime;

print OUTPUT "Attraction model estimation script started at $starttime\n";
print STDERR "Attraction model estimation script started at $starttime\n";

#<INPUT>;  ###Throw out the first line. ###No wait, don't throw out the first line.

my ($sitesmap) = <INPUT>=~/SITES\s+(\S+)/;
<INPUT>;  ##Discard the DEM map, we don't need it for this script, only the cost model calculator.
my ($nullmap) = <INPUT>=~/NULL\s+(\S+)/;

until ((my $record=<INPUT>)=~/ATTRACTMAPS/) ##Skip the lines that contain cost model inputs.
	{next;}

my (%baseline,%area,%proximity,%habcover,%badcover); ##Initialize global hashes associating map name with type.
until ((my $record=<INPUT>)=~/MODELS/) ##Put the chunk in $record until it contains 'MODELS'
	{
	chomp $record;  ##Remove newline from end of chunk.
	my ($maptype,$hashline)=split /\t/,$record; ###Break the lines down into type and commma-delimited hash initialization lines.
	my @hash_line = split ',',$hashline; ###Split effect levels and mapnames into array units.
#		$hash_line=join ',',@hashline;
	$maptype=lc $maptype;  ### Make $maptype all lower case.
#		print STDERR "$maptype\t$hash_line\n";
	for (my $i=0;$i<scalar(@hash_line);$i+=2) ### Go through the array by pairs of basemap class and map file name.
		{
		given ($maptype)
			{   ###given is the new switch-like control structure built into Perl 5.10 and higher.
			when('baseline'){$baseline{lc $hash_line[$i]}=$hash_line[$i+1];}
			when('area'){$area{lc $hash_line[$i]}=$hash_line[$i+1];}
			when('proximity'){$proximity{lc $hash_line[$i]}=$hash_line[$i+1];}
			when('habcover'){$habcover{lc $hash_line[$i]}=$hash_line[$i+1];}
			when('badcover'){$badcover{lc $hash_line[$i]}=$hash_line[$i+1];}
			default {next;}
			} ### given requires 'use v5.10' in the module calls.
#		print STDERR "$hash_line[$i]\t$baseline{$hash_line[$i]}\n";
		}
	}

my $blockline='#################################################';
print OUTPUT "Input map names extracted...\n"
print STDERR "Input map names extracted...\n";
print OUTPUT "$blockline\n";

### This section was copied directly from the cost model calculator.

my $modelheaders = <INPUT>; ### Copy to column headers for the models
chomp $modelheaders;
my @modelcolumns = split /\t/, lc $modelheaders;  #Grab the column titles.
#@modelcolumns = lc @modelcolumns;
my %modelhash;  ##This will store hashes for all models and their components for later use.
shift @modelcolumns; ##Shift out the one for the "Model" column.

while(my $record=<INPUT>)
	{
	chomp $record;
	my ($model,@effects) = split /\t/,$record;  ##Model number will be the modelhash key referencing a hash containing effectname keys and values will be effect levels. Keep the scalar model number and the array of effect levels for each model.
	my $i=0;  ###We're gonna keep track of array elements with this counter.
	foreach(@modelcolumns) #Combine the effect levels into a hash for each model referenced by model number.
		{
		$modelhash{$model}->{"$_"}=lc $effects[$i]; ##Store a hash of model components associated with each model.
#		print STDERR "$_\t$effects[$i]\n";
		$i++;
		}
	}

print OUTPUT "Model characteristics extracted...\n";
print STDERR "Model characteristics extracted...\n";
print "$blockline\n";
###Proximity and Density Effect Calculations

####These calculations are structured by effect with calculations for each of the vehicles to which they pertain.


print OUTPUT "Calculating proximity effects...\n";
print STDERR "Calculating proximity effects...\n";

foreach(keys %proximity) ###I hypothesized proximity effects on attraction only for migratory birds and short dispersing mammals.  These are exponential decay away from specific habitat resources that should add onto generalized exponential decay at some vehicle specific scale to depict gravity wells around croplands.  The idea here is that being close to a resource, regardless of it's area or % habitat cover, makes an area more desirable as a foraging or home range area.
	{
	given (lc $_)
		{
		when ('wetland')
			{
			system "r.mapcalc \"MB_wetland_prox_attract=0-(-0.5*exp(0.003*1200)/(exp(0.003*$proximity{$_})+(exp(0.003*1200))))+0.5\"";
			wait;
			system "r.mapcalc \"MB_wetland_prox_attract=if(isnull(MB_wetland_prox_attract) or MB_wetland_prox_attract<0.1,0.1,MB_wetland_prox_attract)\"";
			wait;
			}
		when ('rodentcover')  ###We already have the density contributing to this.  Instinctively, I think that the proximity and the density probably impact attraction, and proximity probably impacts cost.  Proximity in cost models influences corridor shape and length, but proximity in attraction impacts home range extent around cover.
			{
			system "r.mapcalc \"ST_cover_prox_attract=0-(-0.5*exp(0.05*100)/(exp(0.05*$proximity{$_})+(exp(0.05*100))))+0.5\"";
			wait;
			system "r.mapcalc \"ST_cover_prox_attract=if(isnull(ST_cover_prox_attract) or ST_cover_prox_attract<0.1,0.1,ST_cover_prox_attract)\"";
			wait;
			}
		when('lbroost')
			{
			system "r.mapcalc \"LB_roost_prox_attract=1/(exp(0.008*$proximity{$_}))\""; ##Low...confirmed function shape 13MAR2013
			wait;
			system "r.mapcalc \"LB_roost_prox_attract=if(isnull(LB_roost_prox_attract) or LB_roost_prox_attract<0.1,0.1,LB_roost_prox_attract)\"";
			wait;
			}
		}
	}

print OUTPUT "Proximity calculations done...\n";
print OUTPUT "$blockline\n";
print STDERR "Proximity calculations done...\n";

print OUTPUT "Calculating Habitat Percent Cover Effects...\n";
print STDERR "Calculating Habitat Percent Cover Effects...\n";

foreach (keys %habcover)
	{
	given (lc $_)
		{

		when('lt')   ####Each of these should have an effect calculator for effect x percent habitat coverage.
			{
			system "r.mapcalc \"LT_cover_density_attract=1-(1/exp(0.06*$habcover{$_}))\"";
			wait;
			system "r.mapcalc \"LT_cover_density_attract=if(isnull(LT_cover_density_attract) or LT_cover_density_attract < 0.1,0.1,LT_cover_density_attract)\"";
			wait;
			system "r.mapcalc \"LT_cover_density_half_attract=1-(1/exp(0.03*$habcover{$_}))\"";
			wait;
			system "r.mapcalc \"LT_cover_density_half_attract=if(isnull(LT_cover_density_half_attract) or LT_cover_density_half_attract<0.1,0.1,LT_cover_density_half_attract)\"";
			wait;
			}
		when('st')
			{
			system "r.mapcalc \"ST_cover_density_attract=1.0*(-1+exp(0.2*30)+exp(0.2*$habcover{$_}))/((exp(0.2*$habcover{$_})+1*exp(0.2*30)*10))\"";
			wait;
			system "r.mapcalc \"ST_cover_density_attract=if(isnull(ST_cover_density_attract) or ST_cover_density_attract < 0.1,0.1,ST_cover_density_attract)\"";
			wait;
			system "r.mapcalc \"ST_cover_density_half_attract=1.0*(-1+exp(0.1*30)+exp(0.1*$habcover{$_}))/((exp(0.1*$habcover{$_})+1*exp(0.1*30)*10))\"";
			wait;
			system "r.mapcalc \"ST_cover_density_half_attract=if(isnull(ST_cover_density_half_attract) or ST_cover_density_half_attract<0.1,0.1,ST_cover_density_half_attract)\"";
			wait;
			}
		when('mb')
			{
			system "r.mapcalc \"MB_cover_density_attract=1.0*(-1+exp(0.12*5)+exp(0.12*$habcover{$_}))/((exp(0.12*$habcover{$_})+1*exp(0.12*5)*10))\"";
			wait;
			system "r.mapcalc \"MB_cover_density_attract=if(isnull(MB_cover_density_attract) or MB_cover_density_attract < 0.1,0.1,MB_cover_density_attract)\"";
			wait;
			system "r.mapcalc \"MB_cover_density_half_attract=1.0*(-1+exp(0.06*5)+exp(0.06*$habcover{$_}))/((exp(0.06*$habcover{$_})+1*exp(0.06*5)*10))\"";
			wait;
			system "r.mapcalc \"MB_cover_density_half_attract=if(isnull(MB_cover_density_half_attract) or MB_cover_density_half_attract<0.1,0.1,MB_cover_density_half_attract)\"";
			wait;
			}
		when('lb')
			{
			system "r.mapcalc \"LB_cover_density_attract=if($habcover{$_}>40,1.0,0.025*$habcover{$_})\"";
			wait;
			system "r.mapcalc \"LB_cover_density_attract=if(isnull(LB_cover_density_attract) or LB_cover_density_attract<0.1,0.1,LB_cover_density_attract)\"";
			wait;
			system "r.mapcalc \"LB_cover_density_half_attract=if($habcover{$_}>40,1.0,0.0125*$habcover{$_}\"";
			wait;
			system "r.mapcalc \"LB_cover_density_half_attract=if(isnull(LB_cover_density_half_attract) or LB_cover_density_half_attract<0.1,0.1,LB_cover_density_half_attract\"";
			wait;
			}
		}
	}
print OUTPUT "Habitat cover effect calculations completed...\n";
print OUTPUT "$blockline\n";
print STDERR "Habitat cover effect calculations completed...\n";

print OUTPUT "Calculating bad land cover effects...\n";
print STDERR "Calculating bad land cover effects...\n";
foreach (keys %badcover)
	{
	given (lc $_)
		{	
		when('lt')  ####Each of these cases should have an effect calculator for effect x percent bad habitat coverage.
			{
			system "r.mapcalc \"LT_badcover_attract=0-(-1*exp(0.1*12)/(exp(0.1*$badcover{$_})+(exp(0.1*12)-1)))\"";
			wait;
			system "r.mapcalc \"LT_badcover_attract=if(isnull(LT_badcover_attract) or LT_badcover_attract < 0.01,0.01,LT_badcover_attract)\"";
			wait;
			}  #For LT, > 40% bad habitat cover should be 
		when('st')
			{
			system "r.mapcalc \"ST_badcover_attract=0-(-1.0*(exp(0.15*60))/(exp(0.15*($badcover{$_}))+(exp(0.15*60)-1.0)))\"";
			wait;
			system "r.mapcalc \"ST_badcover_attract=if(isnull(ST_badcover_attract) or ST_badcover_attract < 0.01,0.01,ST_badcover_attract)\"";
			wait;			
			}
		when('mb')
			{
			system "r.mapcalc \"MB_badcover_attract=0-(-1.0*(exp(0.12*55))/((exp(0.12*$badcover{$_})+(exp(0.12*55)-1.0))))\"";
			wait;
			system "r.mapcalc \"MB_badcover_attract=if(isnull(MB_badcover_attract) or MB_badcover_attract < 0.1,0.1,MB_badcover_attract)\"";
			wait;			
			}
		when('lb')
			{
			system "r.mapcalc \"LB_badcover_attract=0-(-1.0*(exp(0.1*15))/((exp(0.1*$badcover{$_}))+(exp(0.1*15)-1.0)))\"";
			wait;
			system "r.mapcalc \"LB_badcover_attract=if(isnull(LB_badcover_attract) or LB_badcover_attract<0.1,0.1,LB_badcover_attract)\"";
			wait;
			}
		}
	}
print OUTPUT "Bad land cover effect calculations complete...\n";
print STDERR "Bad land cover effect calculations complete...\n";
print OUTPUT "$blockline\n";

############################################

####Area effect calculation ################

### This section should be only four r.mapcalc lines: one for each vehicle.  ST and LB may not even need one. 24JUL2012.

print OUTPUT "Calculating area x vehicle effects...\n";

system "r.mapcalc \"LT_area_effect=1*((-1+exp(1))+(1.5*exp($area{'patcharea'})))/(1.5*exp($area{'patcharea'})-1+exp(1)*100)\"";  ####LT
wait;
system "r.mapcalc \"ST_area_effect=if($area{'patcharea'}<=0.04,0.5,1.0)\"";
wait;
system "r.mapcalc \"MB_area_effect=1*((-1+exp(1))+(1.25*exp($area{'patcharea'})))/(1.25*exp($area{'patcharea'})-1+exp(1)*100)\"";  ### Couldn't find any documentation from this, but my own observations are that I haven't seen migratory flocks of Geese in fields smaller than about 2 ha. So, I made it slightly shallower than the LT area effect function. 
wait;
system "r.mapcalc \"LB_area_effect=1.0*((exp($area{'patcharea'}))/(10+exp($area{'patcharea'})))\"";
wait;

print OUTPUT "Calculated area x vehicle effects...\n";
print OUTPUT "$blockline\n";
print STDERR "Calculated area x vehicle effects...\n";

############################################

###Calculate Attraction Models  ###This has contains a lot of miscalls b/c I didn't know the actual relationship between input maps and effect maps when I modified it from the cost map calculator initially.

foreach(keys %modelhash)
	{
	my $modelname = $_;
	my %componenthash=%{$modelhash{$modelname}};
	if($componenthash{'attraction'} eq 'none') ###Logic flow control: if attraction is not 'None' then calculate 'Norm' with the main code.  The other levels for Attraction provide only very specific modifiers to components of 'Norm'
		{
		system "r.mapcalc \"Attract$modelname=$nullmap;\""; 
		wait;
		}
	else
		{
		given($componenthash{'baseline'})
			{
			when('lt')  #Calculate by vehicle class (baseline).  
				{
				given (lc $componenthash{'attraction'})
					{
					when ('norm')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'lt'}*LT_area_effect*LT_badcover_attract*LT_cover_density_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.01,0.01,Attract$modelname)\"";  ###This will bring things back up to the minimum and make sure nulls aren't propagated.  We intentionally set the value very low because some of these areas are places where you would never expect to see these vehicles.
						wait;
						}
					when ('area')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'lt'}*LT_badcover_attract*LT_cover_density_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.01,0.01,Attract$modelname)\"";
						wait;
						}
					when ('hab')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'lt'}*LT_badcover_attract*LT_cover_density_half_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.01,0.01,Attract$modelname)\"";
						wait;
						}
					}
				###Calculate habcover coefficient as half the slope where Attraction='Hab'.  Leave badcover unmodified.
				###Exclude area effect where Attraction='Area' or 'Hab'
				###How do I run proximity away from both internal and external sides of edges?
#				system "r.mapcalc \"$baseline{'lt'}*$area{'lt'}*$habcover{'lt'}*$badcover{'lt'}*$proximity{'lt'}\""; ##Final calculation, basically everything here needs to be replaced with mapnames for effect maps calculated earlier in the script.  The basemaps for all these no longer represent effects, only basemaps from which to calculate effects.
				}
			when('st')
				{
				given (lc $componenthash{'attraction'})
					{
					when ('norm')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'st'}*ST_area_effect*ST_badcover_attract*ST_cover_density_attract*ST_cover_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.01,0.01,Attract$modelname)\"";  ###This will bring things back up to the minimum and make sure nulls aren't propagated.  We intentionally set the value very low because some of these areas are places where you would never expect to see these vehicles.
						wait;
						}
					when ('area')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'st'}*ST_badcover_attract*ST_cover_density_attract*ST_cover_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.01,0.01,Attract$modelname)\"";
						wait;
						}
					when ('hab')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'st'}*ST_badcover_attract*ST_cover_density_half_attract*ST_cover_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.01,0.01,Attract$modelname)\"";
						wait;
						}
					}
				}
			when('mb')
				{
				given (lc $componenthash{'attraction'})
					{
					when ('norm')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'mb'}*MB_area_effect*MB_badcover_attract*MB_cover_density_attract*MB_wetland_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.1,0.1,Attract$modelname)\"";  ###This will bring things back up to the minimum and make sure nulls aren't propagated.  Birds get a higher minimum than terrestrials.
						wait;
						}
					when ('area')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'mb'}*MB_badcover_attract*MB_cover_density_attract*MB_wetland_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.1,0.1,Attract$modelname)\"";
						wait;
						}
					when ('hab')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'mb'}*MB_badcover_attract*MB_cover_density_half_attract*MB_wetland_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.1,0.1,Attract$modelname)\"";
						wait;
						}
					}
				}
			when('lb')
				{
				given (lc $componenthash{'attraction'})
					{
					when ('norm')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'lb'}*LB_area_effect*LB_badcover_attract*LB_cover_density_attract*LB_roost_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.1,0.1,Attract$modelname)\"";  ###This will bring things back up to the minimum and make sure nulls aren't propagated.  Birds get a higher minimum than terrestrials.
						wait;
						}
					when ('area')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'lb'}*LB_badcover_attract*LB_cover_density_attract*LB_roost_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.1,0.1,Attract$modelname)\"";
						wait;
						}
					when ('hab')
						{
						system "r.mapcalc \"Attract$modelname=$baseline{'lt'}*LT_badcover_attract*LT_cover_density_half_attract*LB_roost_prox_attract\"";
						wait;
						system "r.mapcalc \"Attract$modelname=if(isnull(Attract$modelname) or Attract$modelname<0.1,0.1,Attract$modelname)\"";
						wait;
						}
					}

				}
			default{next;print OUTPUT "Attract$modelname not created. Attraction map vehicle undefined.\n";}
			}		
		}
	}

