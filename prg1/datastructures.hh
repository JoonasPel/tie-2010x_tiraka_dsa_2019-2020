// Datastructures.hh

#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <unordered_map>

// Types for IDs
using StopID = long int;
using RegionID = std::string;
using Name = std::string;

// Return values for cases where required thing was not found
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


// This is the class you are supposed to implement

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // Estimate of performance:
    // O(1)
    // Short rationale for estimate:
    // Unordered_mapin .size on vakioaikainen
    int stop_count();

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // Unordered_mapin .clear on lineaarinen
    void clear_all();

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // For-loop O(n), push_back O(1), .first O(1)
    std::vector<StopID> all_stops();

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1).
    bool add_stop(StopID id, Name const& name, Coord xy);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // Sama .find kuin add_stopissa. stops_[id].name pitäisi olla O(1).
    Name get_stop_name(StopID id);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // Sama .find kuin add_stopissa. stops_[id].coord pitäisi olla O(1).
    // Tämä operaatio on lähes identtinen get_stop_name kanssa.
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

private:

    //Tietorakenne pysäkeille.
    struct Stop {
        Name name;
        Coord coord;
        RegionID in_region = NO_REGION;
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

    std::vector<RegionID> find_parents_recursive
    (RegionID id, std::vector<RegionID> vec);

};

#endif // DATASTRUCTURES_HH
