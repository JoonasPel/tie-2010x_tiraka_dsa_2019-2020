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
    // .find on O(n) joskus(huonoin tapaus), keskimäärin Θ(1).
    bool add_stop(StopID id, Name const& name, Coord xy);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). stops_[id].name pitäisi olla O(1)
    Name get_stop_name(StopID id);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). stops_[id].coord pitäisi olla O(1)
    Coord get_stop_coord(StopID id);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // O(n*log n)
    // Short rationale for estimate:
    // Hitain osa on std::sort, joka on O(n*log n)
    std::vector<StopID> stops_alphabetically();

    // Estimate of performance:
    // O(n*log n)
    // Short rationale for estimate:
    // Sama kun ylempi funktio. Hitain osa on std::sort O(n*log n)
    std::vector<StopID> stops_coord_order();

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // For-loop on O(n). Muut toiminnot vakioaikaisia.
    StopID min_coord();

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // For-loop on O(n). Muut toiminnot vakioaikaisia.
    StopID max_coord();

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // for_each on O(n) tehokkuudeltaan.
    std::vector<StopID> find_stops(Name const& name);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1).
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). Nimen vaihto O(1).
    bool change_stop_name(StopID id, Name const& newname);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1).
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). Koordinaatin vaihto O(1).
    bool change_stop_coord(StopID id, Coord newcoord);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). Sama pätee alkion lisäykselle.
    bool add_region(RegionID id, Name const& name);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1)
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). regions_[id].name on O(1)
    Name get_region_name(RegionID id);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // For-loop O(n). Muut toiminnot O(1).
    std::vector<RegionID> all_regions();

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1).
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1).
    // Muut toiminnot O(1) tai keskimäärin Θ(1).
    bool add_stop_to_region(StopID id, RegionID parentid);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // .find on O(n) joskus, keskimäärin Θ(1). Rekursio on O(n).
    // Muut toiminnot ovat O(1) tai keskimäärin Θ(1).
    bool add_subregion_to_region(RegionID id, RegionID parentid);

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // .find on huonoimmillaan O(n), keskimäärin Θ(1). Rekursio on O(n)
    // Käytännössä ohjelmassa ei suuria määriä regioneita lisäillä ainakaan
    // perfteisteissä, joten päästään usein tehokkuuteen O(log n)
    std::vector<RegionID> stop_regions(StopID id);

    // Non-compulsory operations

    // Estimate of performance:
    // Short rationale for estimate:
    void creation_finished();

    // Estimate of performance:
    // O(n)
    // Short rationale for estimate:
    // Loopit ovat lineaarisia ja rekursio on lineaarinen.
    // Huomioidaan, että myös "tuplalooppi" on lineaarinen, koska käytännössä
    // se vain käy kaikki stopit kertaalleen läpi.
    std::pair<Coord, Coord> region_bounding_box(RegionID id);

    // Estimate of performance:
    // O(n*log n)
    // Short rationale for estimate:
    // Loopit ovat lineaarisia, std::sort on (n*log n).
    // Kuitenkin perftesteissä jostain syystä saadaan tehokkuudeksi O(log n).
    // Oletetaan kuitenkin, että oikea tehokkuus on O(n*log n)
    std::vector<StopID> stops_closest_to(StopID id);

    // Estimate of performance:
    // O(n) huonoimmassa tapauksessa, keskimäärin Θ(1).
    // Short rationale for estimate:
    // .find on joskus O(n), keskimäärin Θ(1).
    // Muut toiminnot O(1) tai keskimäärin Θ(1).
    bool remove_stop(StopID id);

    // Estimate of performance:
    // Huonoimmassa tapauksessa O(n^2), keskimäärin Θ(log n)
    // Short rationale for estimate:
    // Tuplalooppi, jossa kahden vektorin kaikkia alkioita verrataan toisiinsa.
    // Huonoimmassa tapauksessa regionit on jaettu kahteen oksaan. Verratut
    // pysäkit löytyy oksien päistä ja yhteinen region on koko puun juuri
    // (tai yhteistä ei ole). Käytännössä tämä on ohjelmassa harvinaista.
    // Keskimääräisesti funktio vaikuttaa olevan logaritminen.
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

};

#endif // DATASTRUCTURES_HH
