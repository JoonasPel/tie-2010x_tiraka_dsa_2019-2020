// Datastructures.hh

#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <unordered_map>
#include <list>
#include <iterator>
#include <stack>
#include <queue>

// Types for IDs
using StopID = long int;
using RegionID = std::string;
using RouteID = std::string;
using Name = std::string;

// Return values for cases where required thing was not found
RouteID const NO_ROUTE = "!!NO_ROUTE!!";
StopID const NO_STOP = -1;
RegionID const NO_REGION = "!!NO_REGION!!";

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();

// Return value for cases where name values were not found
Name const NO_NAME = "!!NO_NAME!!";

// Type for a coordinate (x, y)
struct Coord
{
    int x = NO_VALUE;
    int y = NO_VALUE;
};

// Example: Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); } // Not strictly necessary

// Example: Defining < for Coord so that it can be used
// as key for std::map/set
inline bool operator<(Coord c1, Coord c2)
{
    if (c1.y < c2.y) { return true; }
    else if (c2.y < c1.y) { return false; }
    else { return c1.x < c2.x; }
}

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

// Type for time of day in minutes from midnight (i.e., 60*hours + minutes)
using Time = int;

// Return value for cases where color was not found
Time const NO_TIME = std::numeric_limits<Time>::min();

// Type for a duration of time (in minutes)
using Duration = int;

// Return value for cases where Duration is unknown
Duration const NO_DURATION = NO_VALUE;

// Type for a distance (in metres)
using Distance = int;

// Return value for cases where Duration is unknown
Distance const NO_DISTANCE = NO_VALUE;

//Solmujen värit algoritmejä varten
int const WHITE = 0;
int const GREY = 1;
int const BLACK = 2;

// This is the class you are supposed to implement

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // Estimate of performance:
    // Short rationale for estimate:
    int stop_count();

    // Estimate of performance:
    // Short rationale for estimate:
    void clear_all();

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<StopID> all_stops();

    // Estimate of performance:
    // Short rationale for estimate:
    bool add_stop(StopID id, Name const& name, Coord xy);

    // Estimate of performance:
    // Short rationale for estimate:
    Name get_stop_name(StopID id);

    // Estimate of performance:
    // Short rationale for estimate:
    Coord get_stop_coord(StopID id);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<StopID> stops_alphabetically();

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<StopID> stops_coord_order();

    // Estimate of performance:
    // Short rationale for estimate:
    StopID min_coord();

    // Estimate of performance:
    // Short rationale for estimate:
    StopID max_coord();

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<StopID> find_stops(Name const& name);

    // Estimate of performance:
    // Short rationale for estimate:
    bool change_stop_name(StopID id, Name const& newname);

    // Estimate of performance:
    // Short rationale for estimate:
    bool change_stop_coord(StopID id, Coord newcoord);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // Short rationale for estimate:
    bool add_region(RegionID id, Name const& name);

    // Estimate of performance:
    // Short rationale for estimate:
    Name get_region_name(RegionID id);

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<RegionID> all_regions();

    // Estimate of performance:
    // Short rationale for estimate:
    bool add_stop_to_region(StopID id, RegionID parentid);

    // Estimate of performance:
    // Short rationale for estimate:
    bool add_subregion_to_region(RegionID id, RegionID parentid);

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<RegionID> stop_regions(StopID id);

    // Non-compulsory operations

    // Estimate of performance:
    // Short rationale for estimate:
    void creation_finished();

    // Estimate of performance:
    // Short rationale for estimate:
    std::pair<Coord, Coord> region_bounding_box(RegionID id);

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<StopID> stops_closest_to(StopID id);

    // Estimate of performance:
    // Short rationale for estimate:
    bool remove_stop(StopID id);

    // Estimate of performance:
    // Short rationale for estimate:
    RegionID stops_common_region(StopID id1, StopID id2);

    // Phase 2 operations

    // Estimate of performance:
    // Θ(n)
    // Short rationale for estimate:
    // For-loop suoritetaan N kertaa ( N = reittien määrä)
    // Vektorin push_back on vakioaikainen, paitsi joskus lineaarinen jos
    // tapahtuu "uudelleenvaraamista". (reallocation)
    std::vector<RouteID> all_routes();

    // Estimate of performance:
    // Θ(n)
    // Short rationale for estimate:
    // .find on lineaarinen, loopit lineaarisia
    bool add_route(RouteID id, std::vector<StopID> stops);

    // Estimate of performance:
    // Θ(n)
    // Short rationale for estimate:
    // .find lineaarinen, unordered_mapista haku keskimäärin Θ(1).
    std::vector<std::pair<RouteID, StopID>> routes_from(StopID stopid);

    // Estimate of performance:
    // Θ(n)
    // Short rationale for estimate:
    // .find lineaarinen, loop lineaarinen, push_back vakioaikainen.
    std::vector<StopID> route_stops(RouteID id);

    // Estimate of performance:
    // O(n^2)
    // Short rationale for estimate:
    // .clear mapille sekä vektorille on O(n).
    // Koska toinen .clear on loopissa, saadaan n^2.
    void clear_routes();

    // Estimate of performance:
    // Θ(n). n = pysäkit + reitit
    // Short rationale for estimate:
    // Alun virhetarkastelut ja lopun while-loop lineaarisia.
    // BFS on O(V+E), jossa V= pysäkkien määrä ja E= reittien määrä
    // perftestin perusteellakin lineaarisuus pitää paikkaansa, tosin
    // isoissa alkiomäärissä (kymmeniätuhansia), vakiokerroin alkaa kasvamaan.
    std::vector<std::tuple<StopID, RouteID, Distance>> journey_any(StopID fromstop, StopID tostop);

//    // Non-compulsory operations

    // Estimate of performance:
    // Θ(n). n = pysäkit + reitit
    // Short rationale for estimate:
    //Tämä kutsuu funktiota journey_any
    std::vector<std::tuple<StopID, RouteID, Distance>> journey_least_stops(StopID fromstop, StopID tostop);

    // Estimate of performance:
    // Θ(n). n = pysäkit
    // Short rationale for estimate:
    // Alun virhetarkastelu lineaarinen, lopun loop lineaarinen
    // DFS O(V), jossa V= pysäkkien määrä.
    std::vector<std::tuple<StopID, RouteID, Distance>> journey_with_cycle(StopID fromstop);

    // Estimate of performance:
    // Θ(n). n = pysäkit + reitit
    // Short rationale for estimate:
    // A* algoritmi on O(V+E), V= pysäkkien määrä ja E= reittien määrä
    // Perftest tukee tätä ja tässäkin algoritmissa vakiokerroin tuntuu
    // kasvavan alkiomäärien noustessa kymmeniintuhansiin.
    std::vector<std::tuple<StopID, RouteID, Distance>> journey_shortest_distance(StopID fromstop, StopID tostop);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // .find lineaarinen, push_back vakioaikainen.
    bool add_trip(RouteID routeid, const std::vector<Time> &stop_times);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // .find/.find_if lineaarisia, loop lineaarinen.
    std::vector<std::pair<Time, Duration> > route_times_from(RouteID routeid, StopID stopid);

    // Estimate of performance:
    // Short rationale for estimate:
    std::vector<std::tuple<StopID, RouteID, Time>> journey_earliest_arrival(StopID fromstop, StopID tostop, Time starttime);

    // Estimate of performance:
    // Short rationale for estimate:
    void add_walking_connections(); // Note! This method is completely optional, and not part of any testing

private:

    //Tietorakenne pysäkeille.
    struct Stop {
        Name name;
        Coord coord;
        StopID id; //Lisätty prg2:ssa
        RegionID in_region = NO_REGION;

        //Seuraavat lisätty prg2:ssa:
        std::vector<std::pair<RouteID,StopID>> next_stops = {};
        int color = 0; //"Väri" algoritmeja varten.
        StopID from = NO_STOP; //Algoritmia varten "leivänmuru".
        RouteID with_route = NO_ROUTE; //Mitä reittiä pitkin tultiin.
        Distance dist_from_start = 0; //Matka tästä pysäkistä aloituspysäkkiin.
    };

    std::unordered_map<StopID, Stop> stops_ = {};

    //"Puu" tietorakenne regioneille.
    struct region_node {

        Name name;
        std::vector<RegionID> subregions = {};
        std::vector<StopID> stops = {};
        bool is_subregion = false;
        RegionID parentid = NO_REGION;

    };
    std::unordered_map<RegionID, region_node> regions_;

    double distance_from_origo(Coord coord);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // Suoritetaan niin monta kertaa kuin regioneita löytyy.
    // n on regioneiden määrä samassa oksassa.
    std::vector<RegionID> find_parents_recursive
    (RegionID id, std::vector<RegionID> vec);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // Suoritetaan niin monta kertaa kuin regioneita löytyy.
    // n on regioneiden määrä samassa oksassa.
    std::vector<RegionID> find_children_recursive
    (RegionID id, std::vector<RegionID> vec);


    // Harjoitustyön 2-vaiheen toteutukset

    std::unordered_map<RouteID, std::vector<Stop*>> routes_ = {};

    std::unordered_map<RouteID, std::vector<std::vector<Time>>> trips_ = {};

    //Stop-structiin lisätty tietoja, joita algoritmit tarvitsevat.

    //Määritelty värit solmuille(int). 0 = valkoinen, 1 = harmaa, 2 = musta

    //Laskee kahden pysäkin välisen etäisyyden.
    int stops_distance(Coord coord1, Coord coord2);

    //Nollaa pysäkkien värin, matkan lähtöpisteestä ja "leivänmurun".
    void reset_stops(bool dijkstra);

    //A*-algoritmin käyttöön
    bool relaxA(StopID u, StopID v, StopID tostop,
               std::priority_queue<std::pair<Distance,StopID>>& queue);

};

#endif // DATASTRUCTURES_HH
