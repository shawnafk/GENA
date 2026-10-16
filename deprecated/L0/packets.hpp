#ifndef PACKETS_HPP
#define PACKETS_HPP
#include <utility>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include "structures.hpp"
using namespace std;
uint16_t B2L16(uint16_t num);
uint32_t B2L32(uint32_t num);
uint64_t B2L64(uint64_t num);
int signed_7bits(uint8_t bits);
void seek_ep1(Ep1 & ep, double & HV_Deflect_N, double & HV_Deflect_P, double & HV_Cfilm, double & HV_E, double & HV_MCP1, double & HV_MCP2,\
    unsigned short & Event_Num, unsigned short & Event_Eff, unsigned short & Event_Aux);
void seek_ep2(Ep2 &ep, double & utc_time, double & ant_angle, double & ant_azimuth, double & Star_angular_momentum ,double &  Star_Z_Solar_angle ,int & Shadow ,\
    double & Moon_J2000_q1 ,double & Moon_J2000_q2 ,double & Moon_J2000_q3 ,double & Moon_J2000_Star_X ,double & Moon_J2000_Star_Y ,double & Moon_J2000_Star_Z);
void seek_event(Event& ev, unsigned int& event_D_mark, unsigned int& event_Incid,unsigned int& event_D_incid_m, unsigned int& event_D_incid_z_m, int& event_Times, int& event_Timeus, unsigned int& event_EplHe, int event_fly_time[][2], int event01_data[][2],int event03_data[][2], int event05_data[][2],int event07_data[][2]);
//In packet 
//Out L0
void parsePacket(Packet& packet, GenaL0& l0);
vector<int64_t> get_packet_header_addr(uint8_t * buffer, long long len_range[2]);
vector<int64_t> find_2064(const vector<int64_t>& a, int64_t L);
#endif
