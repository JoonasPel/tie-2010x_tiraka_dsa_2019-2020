// Datastructures.cc

#include "datastructures.hh"

#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <QDebug>
#include <climits>

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
    routes_.clear();
}

std::vector<StopID> Datastructures::all_stops()
{
    std::vector<StopID> stop_IDs =  {};

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

    return stops_[id].name;
}

Coord Datastructures::get_stop_coord(StopID id)
{
    if(stops_.find(id) == stops_.end()) { return NO_COORD; }

    return stops_[id].coord;
}

std::vector<StopID> Datastructures::stops_alphabetically()
{
    std::vector<std::pair<Name,StopID>> stop_pairs;
    std::vector<StopID> ordered_stops = {};

    for (std::pair stop : stops_) {

        stop_pairs.push_back(make_pair(stop.second.name, stop.first));
    }

    std::sort(stop_pairs.begin(), stop_pairs.end(), [](auto &left, auto&right) {
        return left.first < right.first; });

    for(auto stop : stop_pairs) {
        ordered_stops.push_back(stop.second);
    }

    return ordered_stops;
}

std::vector<StopID> Datastructures::stops_coord_order()
{
    //Pitkälti hyödynnetty yllä olevaa stops_alphabetically funktiota

    std::vector<std::pair<Coord,StopID>> stop_pairs;
    std::vector<StopID> ordered_stops = {};

    for (std::pair stop : stops_) {

        stop_pairs.emplace_back(stop.second.coord,stop.first);
    }

    //En käytä tässä funktiota distance_from_origo, koska yksinkertaisempi ilman.
    std::sort(stop_pairs.begin(), stop_pairs.end(), [](auto &left, auto&right) {
        return sqrt(pow(left.first.x,2) + pow(left.first.y,2)) <
                sqrt(pow(right.first.x,2) + pow(right.first.y,2));
    });

    for(auto stop : stop_pairs) {
        ordered_stops.push_back(stop.second);
    }

    return ordered_stops;
}

StopID Datastructures::min_coord()
{
    if(stops_.size() == 0) { return NO_STOP;}

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
}

StopID Datastructures::max_coord()
{
    if(stops_.size() == 0) { return NO_STOP;}

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

    for (std::pair region : regions_) {
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

    //Tarkistetaan ettei synny syklejä.
    std::vector<RegionID> sub_regions = {};
    sub_regions = find_parents_recursive(parentid, sub_regions);
    if(std::find(sub_regions.begin(), sub_regions.end(), id)
            != sub_regions.end()) { return false; }

    //Lisätään subregion ja muut tarpeelliset tiedot.
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
    //Jos regionia ei ole olemassa.
    if(regions_.find(id) == regions_.end()) {return {{NO_COORD},{NO_COORD}};}


    //Kaikki regionit tallennetaan tähän. Järjestyksellä ei ole väliä.
    std::vector<RegionID> vec = {};

    //Tähän haetaan regionin ja sen subregioneiden pysäkkien koordinaatit.
    std::vector<Coord> coords;

    vec = find_children_recursive(id, vec); //Subregionit tallennetaan
    vec.push_back(id); //Pääregion tallennetaan

    //Käydään regionit ja niiden stopit läpi tallentaen koordinaatit.
    for(RegionID reg_id : vec) {

      for(StopID stop_id : regions_[reg_id].stops) {

          coords.push_back(stops_[stop_id].coord);
      }
    }

    //Jos ei löytynyt ollenkaan pysäkkejä(koordinaatteja), lopetetaan suoritus.
    if(coords.size() == 0) {return {{NO_COORD},{NO_COORD}};}

    //Alkuarvot, josta lähdetään vertaamaan muihin.
    int min_x = coords[0].x;
    int max_x = coords[0].x;
    int min_y = coords[0].y;
    int max_y = coords[0].y;

    //Etsitään x- ja y-koordinaattien maksimit ja minimit.
    //Indeksin 0 pysäkki käydään turhaan, mutta ei vaikuta paljoakaan.
    for(Coord coord : coords) {

        if(coord.x < min_x) { min_x = coord.x;}
        if(coord.x > max_x) { max_x = coord.x;}
        if(coord.y < min_y) { min_y = coord.y;}
        if(coord.y > max_y) { max_y = coord.y;}
    }

    return {{min_x, min_y},{max_x, max_y}};
}

std::vector<StopID> Datastructures::stops_closest_to(StopID id)
{
    if(stops_.find(id) == stops_.end()) { return {NO_STOP}; }

    //Tähän tallentuu lähimmät pysäkit.
    std::vector<StopID> closest = {};

    //Koordinaatit, joihin muita pysäkkejä verrataan.
    int x0 = stops_[id].coord.x;
    int y0 = stops_[id].coord.y;


    std::vector<std::pair<Coord,StopID>> stop_pairs;

    for (std::pair stop : stops_) {

        stop_pairs.emplace_back(stop.second.coord,stop.first);
    }

    std::sort(stop_pairs.begin(), stop_pairs.end(), [x0,y0](auto &left, auto&right) {
        return sqrt(pow((left.first.x-x0),2) + pow((left.first.y-y0),2)) <
                sqrt(pow((right.first.x-x0),2) + pow((right.first.y-y0),2));
    });

    //Aloitetaan ykkösestä, koska indeksi 0 on parametrina saatu pysäkki.
    //Lopetetaan, jos pysäkit loppuu kesken tai ollaan saatu 5 lähintä.
    for(unsigned i=1; i < stop_pairs.size() and i < 6; ++i) {

        closest.push_back(stop_pairs[i].second);
    }

    return closest;
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
    //Onko pysäkit olemassa
    if(stops_.find(id1) == stops_.end()) { return {NO_REGION}; }
    if(stops_.find(id2) == stops_.end()) { return {NO_REGION}; }

    //Haetaan pysäkkien regionit käyttäen valmista toteutusta.
    std::vector<RegionID> stop1_regions = stop_regions(id1);
    std::vector<RegionID> stop2_regions = stop_regions(id2);

    for(RegionID region_1  : stop1_regions) {

        for(RegionID region_2 : stop2_regions) {

            if(region_1 == region_2) {
                return region_1;
            }
        }
    }

    //Ei löytynyt samaa
    return {NO_REGION};
}

std::vector<RouteID> Datastructures::all_routes()
{
    if(routes_.empty()) { return {NO_ROUTE}; }

    std::vector<RouteID> routeids = {};
    for(std::pair routeid : routes_) {
        routeids.push_back(routeid.first);
    }
    return routeids;
}

bool Datastructures::add_route(RouteID id, std::vector<StopID> stops)
{   
    //Virhetarkastelut puuttuu

    routes_[id] = stops;

    //Lisätään myös pysäkeille tieto, että reitti kulkee sen kautta.
    //Ei lisätä päätepysäkille, koska reitti ei jatku siitä. (siksi size-1)
    for(unsigned int i=0 ; i< stops.size()-1; i++) {
        stops_[stops[i]].stop_routes.push_back(id);

        //Pysäkeille tiedot, mihin pysäkkeihin niistä pääsee.
        StopID from = stops[i];
        StopID to = stops[i+1];
        stops_[from].next_stops.push_back({id,to});
    }

    return true;
}

std::vector<std::pair<RouteID, StopID>> Datastructures::routes_from(StopID stopid)
{
    std::vector<std::pair<RouteID, StopID>> route_stop = {};

    for(RouteID route : stops_[stopid].stop_routes) {

        auto it = std::find(routes_[route].begin(),
                            routes_[route].end(), stopid);
        auto index = std::distance(routes_[route].begin(), it);

        StopID next_stop = routes_[route].at(index+1);
        route_stop.push_back({route,next_stop});

    }

    return route_stop;


    //return {{NO_ROUTE, NO_STOP}};
}

std::vector<StopID> Datastructures::route_stops(RouteID id)
{
    //Jos reittiä ei löydy.
    if(routes_.find(id) == routes_.end()) { return {NO_STOP}; }

    return routes_[id];
}

void Datastructures::clear_routes()
{
    routes_.clear();
}

std::vector<std::tuple<StopID, RouteID, Distance>> Datastructures::journey_any(StopID fromstop, StopID tostop)
{
    //Toteutettu BFS-algoritmilla.

    //Tarkistetaan, että pysäkit löytyy.
    if(stops_.find(fromstop) == stops_.end()) {
        return {{NO_STOP, NO_ROUTE, NO_DISTANCE}};
    } else if (stops_.find(tostop) == stops_.end()) {
        return {{NO_STOP, NO_ROUTE, NO_DISTANCE}};
    }

    std::vector<std::tuple<StopID, RouteID, Distance>> final_route = {};

    reset_stops(false); //falsella matkat asetetaan nollaksi.

    std::vector<StopID> BFS_queue = {};
    bool found_tostop = false;

    stops_[fromstop].color = 1;
    BFS_queue.push_back(fromstop);

    while(!BFS_queue.empty()) {

        StopID id = BFS_queue[0];
        auto *start_ptr = &stops_[id];

        for(auto neighbour : start_ptr->next_stops) {

            auto *neighb_ptr = &stops_[neighbour.second];

            if(neighb_ptr->color == 0) {
                neighb_ptr->color = 1;

                neighb_ptr->dist_from_start = start_ptr->dist_from_start +
                        stops_distance(start_ptr->coord,neighb_ptr->coord);

                neighb_ptr->from = id;
                neighb_ptr->with_route = neighbour.first;

                //Tarkistaa löydettiinkö etsitty määränpää (tostop)
                if(neighbour.second == tostop) {
                    found_tostop = true;
                    break;
                }
                BFS_queue.push_back(neighbour.second);
            }
        }
        if(found_tostop) { break; }

        start_ptr->color = 2;
        BFS_queue.erase(BFS_queue.begin());
    }
    //Jos tostopia ei löytynyt (ei ollut reittiä) se on valkoinen.
    if(stops_[tostop].color == 0) {return {};}

    StopID temp = tostop;
    RouteID route = NO_ROUTE;
    auto *stop_ptr = &stops_[temp];

    while(temp != fromstop) {
        final_route.push_back({temp, route, stop_ptr->dist_from_start});
        route = stop_ptr->with_route;
        temp = stop_ptr->from;
        stop_ptr = &stops_[temp];
    }
    final_route.push_back({temp, route, stop_ptr->dist_from_start});
    std::reverse(final_route.begin(),final_route.end());

    return final_route;
}

std::vector<std::tuple<StopID, RouteID, Distance>> Datastructures::journey_least_stops(StopID fromstop, StopID tostop)
{
    //journey_any funktiossa käytin BFS:sää, joka sopii tähänkin hyvin.
    //Kutsutaan sitä vain tässä, koska turha copypastea.
    return journey_any(fromstop,tostop);
}

std::vector<std::tuple<StopID, RouteID, Distance>> Datastructures::journey_with_cycle(StopID fromstop)
{   
    //Toteutettu DFS-algoritmilla.

    //Tarkistetaan, että pysäkki löytyy.
    if(stops_.find(fromstop) == stops_.end()) {
        return {{NO_STOP, NO_ROUTE, NO_DISTANCE}};
    }

    std::vector<std::tuple<StopID, RouteID, Distance>> final_route = {};
    bool cycle_found = false;
    StopID last_before_cycle = NO_STOP; //vika pysäkki ennen sykliä
    RouteID route_to_cycle = NO_ROUTE;
    StopID stop_2nd_time = NO_STOP; //pysäkki joka aiheuttaa syklin.
    Distance total_dist = 0;

    reset_stops(false); //falsella matkat asetetaan nollaksi.

    std::stack<StopID> DFS_stack = {};

    DFS_stack.push(fromstop);

    while(!DFS_stack.empty()) {

        if(cycle_found) {break;}

        StopID id = DFS_stack.top();
        auto *start_ptr = &stops_[id];
        DFS_stack.pop();

        if(start_ptr->color == 0) {
            start_ptr->color = 1;
            DFS_stack.push(id);

            for(auto neighbour : start_ptr->next_stops) {

                auto *neighb_ptr = &stops_[neighbour.second];

                if(neighb_ptr->color == 0) {

                    neighb_ptr->dist_from_start = start_ptr->dist_from_start +
                            stops_distance(start_ptr->coord,neighb_ptr->coord);
                    neighb_ptr->from = id;
                    neighb_ptr->with_route = neighbour.first;

                    DFS_stack.push(neighbour.second);

                } else if (neighb_ptr->color == 1) {
                    //Otetaan ylös syklin aiheuttaman pysäkin tiedot, koska
                    //ei voida tallentaa pysäkkiin suoraan.
                    //(Ylikirjoittaisi tiedot ensimmäisestä kerrasta.)
                    last_before_cycle = id;
                    route_to_cycle = neighbour.first;
                    stop_2nd_time = neighbour.second;
                    total_dist = start_ptr->dist_from_start +
                            stops_distance(start_ptr->coord,neighb_ptr->coord);

                    cycle_found = true;
                    break;
                }
            }
        } else {
            start_ptr->color = 2;
        }
    }

    if(!cycle_found) {return {};}

    StopID temp = last_before_cycle;
    RouteID route = route_to_cycle;
    auto *stop_ptr = &stops_[temp];

    final_route.push_back({stop_2nd_time, NO_ROUTE, total_dist});
    while(temp != fromstop) {
        final_route.push_back({temp, route, stop_ptr->dist_from_start});
        route = stop_ptr->with_route;
        temp = stop_ptr->from;
        stop_ptr = &stops_[temp];
    }
    final_route.push_back({temp, route, stop_ptr->dist_from_start});
    std::reverse(final_route.begin(),final_route.end());

    return final_route;
}

std::vector<std::tuple<StopID, RouteID, Distance>> Datastructures::journey_shortest_distance(StopID fromstop, StopID tostop)
{
    //Tarkistetaan, että pysäkit löytyy.
    if(stops_.find(fromstop) == stops_.end()) {
        return {{NO_STOP, NO_ROUTE, NO_DISTANCE}};
    } else if (stops_.find(tostop) == stops_.end()) {
        return {{NO_STOP, NO_ROUTE, NO_DISTANCE}};
    }

    reset_stops(true); //truella matkat asetetaan "äärettömäksi".

    std::vector<std::tuple<StopID, RouteID, Distance>> final_route = {};
    std::priority_queue<std::pair<Distance,StopID>> prio_queue = {};

    auto *from_ptr = &stops_[fromstop];
    from_ptr->color = 1;
    from_ptr->dist_from_start = 0;
    prio_queue.push({from_ptr->dist_from_start,fromstop});
    bool stop = false;

    while(!prio_queue.empty()) {

        //Tällä loopilla poistetaan dublikaatit/valitaan uusi tila.
        //Tarkistaa myös tuleeko jonosta ulos tostop. Jos tulee, algo valmis.
        while(true) {

            std::pair<Distance,StopID> temp = prio_queue.top();
            temp.first = -temp.first; //Muutetaan negatiivisesta positiiviseksi
            if(temp.first > stops_[temp.second].dist_from_start){
                //Löytyi dublikaatti pysäkille, jolle on jo lyhyempi reitti.
                prio_queue.pop();
            } else {
                if(temp.second == tostop) {stop = true;}
                break;
            }
            if(prio_queue.empty()) {
                stop = true;
                break;
            }
        }

        //Jos priority_queue tyhjeni poistettaessa dublikaatteja.
        //Tai jos löydettiin maali eli tostop.
        if(stop) {break;}

        StopID id = prio_queue.top().second;
        auto *start_ptr = &stops_[id];
        prio_queue.pop();

        for(auto neighbour : start_ptr->next_stops) {

            auto *neighb_ptr = &stops_[neighbour.second];

            if(neighb_ptr->color == 0) {
                neighb_ptr->color = 1;
                //Lisätään arvo negatiivisena, koska jono antaa suurimman ensin.
                prio_queue.push({-neighb_ptr->dist_from_start,neighbour.second});
            }
            if(relax(id, neighbour.second,prio_queue)) {
                //Löytyi lyhyempi reitti, päivitetään myös route oikeaksi.
                neighb_ptr->with_route = neighbour.first;
            }          
    }
        start_ptr->color = 2;
    }
    //Jos tostopia ei löytynyt (ei ollut reittiä) se on valkoinen.
    if(stops_[tostop].color == 0) {return {};}

    StopID temp = tostop;
    RouteID route = NO_ROUTE;
    auto *stop_ptr = &stops_[temp];

    final_route.push_back({temp, NO_ROUTE, stop_ptr->dist_from_start});
    while(temp != fromstop) {
        route = stop_ptr->with_route;
        temp = stop_ptr->from;
        stop_ptr = &stops_[temp];
        final_route.push_back({temp, route, stop_ptr->dist_from_start});
    }
    std::reverse(final_route.begin(),final_route.end());

    return final_route;
}

bool Datastructures::add_trip(RouteID routeid, std::vector<Time> const& stop_times)
{
    // Replace this comment and the line below with your implementation
    return false;
}

std::vector<std::pair<Time, Duration>> Datastructures::route_times_from(RouteID routeid, StopID stopid)
{
    // Replace this comment and the line below with your implementation
    return {{NO_TIME, NO_DURATION}};
}

std::vector<std::tuple<StopID, RouteID, Time> > Datastructures::journey_earliest_arrival(StopID fromstop, StopID tostop, Time starttime)
{
    // Replace this comment and the line below with your implementation
    return {{NO_STOP, NO_ROUTE, NO_TIME}};
}

void Datastructures::add_walking_connections()
{
    // Replace this comment and the line below with your implementation
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

std::vector<RegionID> Datastructures::find_children_recursive(RegionID id, std::vector<RegionID> vec)
{
    //Jos regionilla ei ole subregioneita
    if(regions_[id].subregions.size() == 0) {
        return vec;
    } else {

        for(RegionID subregion : regions_[id].subregions) {
            vec.push_back(subregion);
            find_children_recursive(subregion, vec);
        }
    }
    return vec;
}

int Datastructures::stops_distance(Coord coord1, Coord coord2)
{
    unsigned int distance =
            sqrt(pow((coord1.x-coord2.x),2) + pow((coord1.y-coord2.y),2));
    return distance;
}

void Datastructures::reset_stops(bool dijkstra)
{
    Distance value = 0;

    //Dijkstra vaatii matkat "äärettömiksi" eikä nollaksi kuten esim. BFS/DFS
    if(dijkstra) { value = INT_MAX;}

    for(auto &stop : stops_) {
        stop.second.color = 0;
        stop.second.dist_from_start = value;
        stop.second.from = NO_STOP;
    }
}

bool Datastructures::relax(StopID u, StopID v,
                           std::priority_queue<std::pair<Distance,StopID>>& queue)
{
    //Osoittimet vähentämään indeksoinnin määrää.
    auto *ptr_v = &stops_[v];
    auto *ptr_u = &stops_[u];
    Distance u_v_distance = stops_distance(ptr_v->coord,ptr_u->coord);

    if(ptr_v->dist_from_start > ptr_u->dist_from_start + u_v_distance) {
        ptr_v->dist_from_start = ptr_u->dist_from_start + u_v_distance;
        ptr_v->from = u;

        //Lisätään priority_queueen pysäkki v uudella isommalla prioriteetilla.
        queue.push({-ptr_v->dist_from_start, v});

        return true; //Löydettiin lyhyempi reitti
    }
    return false;
}

