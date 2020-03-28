// Datastructures.cc

#include "datastructures.hh"

#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

// Modify the code below to implement the functionality of the class.
// Also remove comments from the parameter names when you implement
// an operation (Commenting out parameter name prevents compiler from
// warning about unused parameters on operations you haven't yet implemented.)

Datastructures::Datastructures()
{
    // Replace this comment with your implementation
}

Datastructures::~Datastructures()
{
    // Replace this comment with your implementation
}

int Datastructures::stop_count()
{
    return stops_.size();
}

void Datastructures::clear_all()
{
    stops_.clear();
    regions_.clear();
}

std::vector<StopID> Datastructures::all_stops()
{
    std::vector<StopID> stop_IDs;

    for (std::pair stop : stops_) {
        stop_IDs.push_back(stop.first);
    }

    return stop_IDs;
}

bool Datastructures::add_stop(StopID id, const Name& name, Coord xy)
{
    //Tarkistetaan onko samalla ID:llä jo pysäkki.
    if(stops_.find(id) != stops_.end()) { return false; }

    stops_[id] = {name,xy};
    return true;
}

Name Datastructures::get_stop_name(StopID id)
{
    if(stops_.find(id) == stops_.end()) { return NO_NAME; }

    std::string stop_name =  stops_[id].name;
    return stop_name;
}

Coord Datastructures::get_stop_coord(StopID id)
{
    if(stops_.find(id) == stops_.end()) { return NO_COORD; }

    Coord stop_coord = stops_[id].coord;
    return stop_coord;
}

std::vector<StopID> Datastructures::stops_alphabetically()
{
    std::vector<std::pair<Name,StopID>> stop_pairs;
    std::vector<StopID> alphabetical_stops;

    for (std::pair stop : stops_) {

        stop_pairs.push_back(make_pair(stop.second.name, stop.first));
    }

    std::sort(stop_pairs.begin(), stop_pairs.end(), [](auto &left, auto&right) {
        return left.first < right.first; });

    for(auto stop : stop_pairs) {
        alphabetical_stops.push_back(stop.second);
    }

    return alphabetical_stops;
    //return {NO_STOP};
}

std::vector<StopID> Datastructures::stops_coord_order()
{
    //Pitkälti hyödynnetty yllä olevaa stops_alphabetically funktiota

    std::vector<std::pair<Coord,StopID>> stop_pairs;
    std::vector<StopID> alphabetical_stops;

    for (std::pair stop : stops_) {

        stop_pairs.emplace_back(stop.second.coord,stop.first);
    }

    //En käytä tässä funktiota distance_from_origo, koska yksinkertaisempi ilman.
    std::sort(stop_pairs.begin(), stop_pairs.end(), [](auto &left, auto&right) {
        return sqrt(pow(left.first.x,2) + pow(left.first.y,2)) <
                sqrt(pow(right.first.x,2) + pow(right.first.y,2));
    });

    for(auto stop : stop_pairs) {
        alphabetical_stops.push_back(stop.second);
    }

    return alphabetical_stops;

    //return {NO_STOP};
}

StopID Datastructures::min_coord()
{
    StopID min_id;
    double min_distance;
    double temp_distance;

    /*
     stops_.begin() ottaa satunnaisen stopin, jonka avulla
     voidaan alkaa verrata stoppeja ja etsiä pienin coord.
    */
    auto random_stop = stops_.begin();
    min_id = random_stop->first;
    min_distance = distance_from_origo(random_stop->second.coord);

    for(std::pair stop : stops_) {

        temp_distance = distance_from_origo(stop.second.coord);

        if(temp_distance < min_distance) {

            min_id = stop.first;
            min_distance = temp_distance;
        }
    }

    return min_id;
    //return NO_STOP;
}

StopID Datastructures::max_coord()
{
    StopID max_id;
    double max_distance;
    double temp_distance;

    auto random_stop = stops_.begin();
    max_id = random_stop->first;
    max_distance = distance_from_origo(random_stop->second.coord);

    for(std::pair stop : stops_) {

        temp_distance = distance_from_origo(stop.second.coord);

        if(temp_distance > max_distance) {

            max_id = stop.first;
            max_distance = temp_distance;
        }
    }

    return max_id;
    //return NO_STOP;
}

std::vector<StopID> Datastructures::find_stops(Name const& name)
{
    std::vector<StopID> found_stops = {};

    std::for_each(stops_.begin(), stops_.end(), [&name, &found_stops](auto &it){
        if (it.second.name == name) {
            found_stops.push_back(it.first);
        }});

    return found_stops;
}

bool Datastructures::change_stop_name(StopID id, const Name& newname)
{
    if(stops_.find(id) == stops_.end()) { return false; }

    stops_[id].name = newname;

    return true;
}

bool Datastructures::change_stop_coord(StopID id, Coord newcoord)
{
    if(stops_.find(id) == stops_.end()) { return false; }

    stops_[id].coord = newcoord;
    return true;

}

bool Datastructures::add_region(RegionID id, const Name &name)
{
    if(regions_.find(id) != regions_.end()) { return false; }

    regions_[id] = {name};
    return true;
}

Name Datastructures::get_region_name(RegionID id)
{
    if(regions_.find(id) == regions_.end()) { return NO_NAME; }

    return regions_[id].name;
}

std::vector<RegionID> Datastructures::all_regions()
{
    std::vector<RegionID> region_IDs = {};

    for (auto region : regions_) {
        region_IDs.push_back(region.first);
    }
    return region_IDs;
}

bool Datastructures::add_stop_to_region(StopID id, RegionID parentid)
{
    //Virhetarkastelut
    if(stops_.find(id) == stops_.end()) { return false; }
    if(regions_.find(parentid) == regions_.end()) { return false; }
    if(stops_[id].in_region != NO_REGION) { return false; }

    regions_[parentid].stops.push_back(id);
    stops_[id].in_region = parentid;
    return true;
}

bool Datastructures::add_subregion_to_region(RegionID id, RegionID parentid)
{
    //Virhetarkastelut
    if(regions_.find(id) == regions_.end()) { return false; }
    if(regions_.find(parentid) == regions_.end()) { return false; }
    if(regions_[id].is_subregion == true) { return false; }

    regions_[id].is_subregion = true;
    regions_[id].parentid = parentid;

    regions_[parentid].subregions.push_back(id);

    return true;
}

std::vector<RegionID> Datastructures::stop_regions(StopID id)
{
    if(stops_.find(id) == stops_.end()) { return {NO_REGION}; }
    if(stops_[id].in_region == NO_REGION) { return {}; }

    std::vector<RegionID> regions = {};

    RegionID stop_in_region = stops_[id].in_region;
    regions.push_back(stop_in_region);

    regions = find_parents_recursive(stop_in_region,regions);

    return regions;
}

void Datastructures::creation_finished()
{
    // Replace this comment with your implementation
    // You don't have to use this method for anything, if you don't need it
}

std::pair<Coord,Coord> Datastructures::region_bounding_box(RegionID id)
{
    // Replace this comment and the line below with your implementation
    return {NO_COORD, NO_COORD};
}

std::vector<StopID> Datastructures::stops_closest_to(StopID id)
{
    // Replace this comment and the line below with your implementation
    return {NO_STOP};
}

bool Datastructures::remove_stop(StopID id)
{
    //Onko pysäkkiä olemassa.
    if(stops_.find(id) == stops_.end()) { return false; }

    //Region, jossa pysäkki on.
    RegionID parentid = stops_[id].in_region;

    //Jos pysäkki ei ole NO_REGIONISSA.
    if(parentid != NO_REGION) {

        //Pysäkin parent_regionissa oleva vektori, jossa regionin kaikki stopit.
        std::vector<StopID>& stops = regions_[parentid].stops;

        stops.erase(std::remove(stops.begin(), stops.end(), id), stops.end());
    }
    stops_.erase(id);

    return true;
}

RegionID Datastructures::stops_common_region(StopID id1, StopID id2)
{
    // Replace this comment and the line below with your implementation
    return NO_REGION;
}

double Datastructures::distance_from_origo(Coord coord)
{
    return sqrt(pow(coord.x,2) + pow(coord.y,2));
}

std::vector<RegionID> Datastructures::find_parents_recursive(RegionID id, std::vector<RegionID> vec)
{
    //Jos regionilla ei ole parentia.
    if(regions_[id].is_subregion == false) {
        return vec;
    } else {

        RegionID parentid = regions_[id].parentid;
        vec.push_back(parentid);
        return find_parents_recursive(parentid,vec);
    }

}


