#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP
#include <cstdint>
#include <vector>
#define SPEED_C 299792458L

//Instrument parameters
#define X_LONG  150.
#define Y_LONG  60.
#define X_MCP  145.
#define Y_MCP  55.
#define X_CF  90.
#define Y_CF  30.



#define NUM_EVENT_PACKETS 62
struct Header {
    //4
    uint8_t bytes[4];
};

struct __attribute__((packed)) Ep1 {
    //12
    //Two's complement
    uint16_t HVs[6];
    //8
    uint8_t monitors[8];
    //3
    uint8_t temperature[3];
    //1
    uint8_t stat;
    //6
    uint16_t event_count[3];
    //2
    uint8_t instruction_counts[2];
    //total 32 
};

struct __attribute__((packed)) Ep2 {
    //1
    uint8_t version;
    //1
    uint8_t err;
    //1
    //used to decide reuse frame
    uint8_t stat;
    //24
    uint8_t flex[24];
    //5
    uint8_t args[5];
    //total 32
};

struct __attribute__((packed)) Event {
    //1
    uint8_t flag;
    //2
    uint16_t time_ms;
    //2
    uint16_t time_us;
    //3
    uint8_t se_amp[3];
    
    //fly time and pulse time; 
    //each time take 2 Bytes space.
    //tof0 has four groups, identical
    //the first (high) 3 bits denotes channel number
    //the last (low) 13 bits denotes tdc value
    uint16_t c0[4];
    //4
    //tof 1 3 5 6 has 2 groups, marks start (second electron) and end (ENA)
    uint16_t c1[2];
    //4
    uint16_t c3[2];
    //4
    uint16_t c5[2];
    //4
    uint16_t c7[2];
    //total 32
};

struct __attribute__((packed)) Packet {
    //4
    Header header;
    //4 2 + 2
    //seq[0] store group mark ...
    //seq[1] packet data domain length
    uint16_t seq[2];
    //4
    uint32_t time_s;
    //2
    uint16_t time_ms;
    //32
    Ep1 ep1;
    //32
    Ep2 ep2;
    //NUM_EVENT_PACKETS*32
    Event evts[NUM_EVENT_PACKETS];
    //2
    std::uint16_t verification_sum;
    //total 14+32*64+2 
};

struct EventEps {
    double HV_Deflect_N, HV_Deflect_P, HV_Cfilm,HV_E,HV_MCP1,HV_MCP2;
    unsigned short Event_Num, Event_Eff, Event_Aux;
    double utc_time,  ant_angle, ant_azimuth, Star_angular_momentum, Star_Z_Solar_angle,\
            Moon_J2000_q1, Moon_J2000_q2, Moon_J2000_q3;
    int Moon_J2000_Star_X, Moon_J2000_Star_Y, Moon_J2000_Star_Z;
    unsigned int Shadow;
};

struct EventL0 {
    unsigned int prob_flag;
    unsigned int eff_flag;
    unsigned int ABmod;
    unsigned int submod;
    double internal_time;
    unsigned int tof[4];
    unsigned int c1[2], c3[2], c5[2], c7[2];
    unsigned int ep;
    unsigned int isfill;
};

//One structure contains 62 events
struct GenaL0 {
    unsigned int BX, BM, Tim_s, Tim_ms;
    EventEps eps;
    EventL0 evts[NUM_EVENT_PACKETS];
};

struct GenaParameters {
    double Timestamp;
    unsigned int MCP1, MCP2, AX, AY;
    int BX1, BY1, BX2, BY2, BG;
    double HV_Deflect_N, HV_Deflect_P, HV_Cfilm,HV_E,HV_MCP1,HV_MCP2;
    double TA, TB;
    unsigned short Event_Num, Event_Eff,Event_Aux;
    double Star_Z_Solar_angle;
    double Moon_J2000_q1;
    double Moon_J2000_q2;
    double Moon_J2000_q3;
    int Moon_J2000_Star_X;
    int Moon_J2000_Star_Y;
    int Moon_J2000_Star_Z;
};
struct GenaParameters_A {
    double Timestamp;
    unsigned int A0, A1, A2, A3;
};
struct GenaParameters_B {
    double Timestamp;
    unsigned int B, B0, B1, B2, B3;
};



//L1A
//non zero
//time to displacement
//A & B label
//One structure only contain one event
struct GenaL1{
    unsigned int flag;
    unsigned int Tim_s, Tim_ms;
    //engineering parameter
    //use l0 directly
    EventEps eps;
    //tof_0, tof_1, tof_2, tof_3, 
    double fly_time[4];
    double tx1[2], tx2[2], ty1[2], ty2[2];
    double ep_mv;
};

//L2
struct GenaL2{
    //quant code
    unsigned int quant;
    //0: good
    //1: exceed mcp
    //2: ns overflow
    //3: not matched
    //AB flag
    unsigned int flag;
    unsigned int Tim_s, Tim_ms;
    //engineering parameter
    double Star_Z_Solar_angle;
    double Moon_J2000_q1;
    double Moon_J2000_q2;
    double Moon_J2000_q3;
 
    int Moon_J2000_Star_X;
    int Moon_J2000_Star_Y;
    int Moon_J2000_Star_Z;

    //use l0 directly
    //x1 ,y1 ,x2 ,y2
    double fly_time;
    double velocity;
    //0: H, 1: O
    int spec;
    double mass;
    //start coordinate
    double x1;
    double y1;
    //stop coordinate
    double x2;
    double y2;
    //energy
    double energy;
    double ep_mv;
};
template<typename T>
struct data_ary {
    T* data;
    unsigned long long count;
};

typedef std::vector<std::pair<int, int>> vec_2int;
typedef std::pair<vec_2int, vec_2int> idx_pair;
typedef std::vector<int> vidx;
typedef std::pair<vidx,vidx> vidx2;
#endif
