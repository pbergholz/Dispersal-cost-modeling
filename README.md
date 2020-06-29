# Dispersal-cost-modeling
Scripts and data for dispersal cost modeling of E. coli using GRASS GIS.
These scripts were developed to automate raster map calculations in GRASS GIS. It consists of three scripts intended to be run in order.


# Software requirements
Perl v 5.10 or greater
GRASS GIS v 6 or greater.


# Summary
The scripts use the GRASS GIS commands r.mapcalc and r.walk to generate landscape resistance maps (resistance against movement) to reflect the movement preferences of wild animals that may carry E. coli from various source locations to produce fields as sink locations.  The GRASS GIS command r.mapcalc is used to iteratively calculate resistance maps based on user formulated basemaps.  Subsequently, r.walk is used to generate cost-distance values among a set of landscape sites.  These values constitute simulations of the dispersal of E. coli with wild animals under various user-defined models of wildlife movement preferences that can then be tested against real population genetic data from E. coli isolated from farm fields and potential source sites.

# Cost_model_calculator.v3.pl

The resistance/cost models are calculated by summing or multiplying among XX different factors affecting landscape resistance:
1)	the baseline landcover preferences of modeled wildlife (e.g. terrestrial ungulates prefer to move through forest, wetland and scrubland and avoid moving through areas with human activity and/or lack of tall vegetation),
2)	resistance from slope for large terrestrial animals and low flying birds,
3)	the effects of proximity from the modeled animals to the nearest landcover (e.g. penetration by ungulates into areas of short/no vegetation cover from preferred areas like forests),
4)	the effects of density of different landcover types (e.g. ungulates might avoid areas/patches of forest that would otherwise be good movement pathways when those patches wind through a dense urban area), 
5)	the effects of barriers (wide streams and roads) that can block the movement of animals, and
6)	the effects of riparian corridors, which a frequently the only means of linking farms in different parts of the landscape with corridors of tall vegetation.
The cost models, of course, involve some simplifying assumptions.  The scripts model four different types of wild animals, but the cost calculations assume that only one type of modeled wild animal is carrying E. coli among locations in any given resistance surface, and that density of landcover affects both attraction to various locations (see below) and cost of moving through various locations.


# attraction_calculator.pl

A second component of E. coli dispersal with modeled wild animals is the attraction of various areas/patches of the wild animal vehicles.  Attractiveness of a patch affects both the draw of such an area for modeled animals and the residence time of those animals in such an area before they move on to a new part of the landscape.  Both components will affect the deposition and acquisition of E. coli by wild animals as they move across the land.  Attraction of locations affects all E. coli vehicles, but it exerts a more important predictive effect on animals that exhibit low cost of movement from location to location under most/all conditions, e.g. birds, humans or other animals that combine the ability to perceive the landscape at larger scales and essentially unlimited dispersal capabilities.
Attraction is modeled as a gravity-like function with attraction scores for destinations (in this study, cropland for growing fresh produce commodities) and exponential decay of the attraction score over distance.  The attraction score for a cropland area in this modeling script is determined by:
1)	baseline attractiveness of the cropland,
2)	patch area (size) with larger patches of attractive land having higher scores,
3)	percent cover of desirable habitat for wild animals that surrounds the target location, and
4)	percent cover of repulsive landcover types surrounding the target location
