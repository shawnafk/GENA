#include <cstdint>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "structures.hpp"
#include "filter.hpp"
#include "packets.hpp"
#include "extern.hpp"
#include <unistd.h> // 包含 Unix/Linux 头文件
#include <algorithm>

using namespace std;

// from L0 to L2 all filter steps
// functions on L0
// TDC to ns
double tdc2tof(int tdc)
{
	return (tdc <= 94) ? (0) : ((tdc - 94) * 0.081);
	//return (tdc - 94) * 0.081;
}

// Ep to mV
double eplhe2ep(unsigned int eplhe)
{
	return (eplhe == 4095 || eplhe == 0) ? (0) : ((2044 - eplhe) * 100.0 / 123.);
}

// given start packets, and num patch
data_ary<GenaParameters> para_from_raw(Packet *packet, unsigned long long num_patch)
{
	GenaParameters * para = (GenaParameters *)malloc(num_patch * sizeof(GenaParameters));
	auto fill_MCP1 = 0;
	auto fill_MCP2 = 0;
	auto fill_AX = 0;
	auto fill_AY = 0;
	auto fill_BX1 = 0;
	auto fill_BX2 = 0;
	auto fill_BY1 = 0;
	auto fill_BY2 = 0;
	auto fill_BG = 0;

	double fill_Moon_J2000_q1 = -1;
	double fill_Moon_J2000_q2 = -1;
	double fill_Moon_J2000_q3 = -1;
	int fill_Moon_J2000_Star_X = -1;
	int fill_Moon_J2000_Star_Y = -1;
	int fill_Moon_J2000_Star_Z = -1;

	auto fill_Star_Z_Solar_angle = -1;
	
	uint32_t num32;
	uint8_t num8;
	unsigned int A;
	for (long i = 0; i < num_patch; i++)
	{
		unsigned int time_s = B2L32(packet[i].time_s);
		unsigned int time_ms = B2L16(packet[i].time_ms);
		double timestamp = time_s + time_ms/1000.;
		para[i].Timestamp = timestamp; 
		A = packet[i].ep1.temperature[0];
        para[i].TA =   (4100*298.15/(298.15*log((A*0.078125*10000)/(5000* (4.97-A*0.078125)))+4100))-273.15;
		A = packet[i].ep1.temperature[1];
        para[i].TB =   (4100*298.15/(298.15*log((A*0.078125*10000)/(5000* (4.97-A*0.078125)))+4100))-273.15;
		seek_ep1(packet[i].ep1,\
			para[i].HV_Deflect_N,\
			para[i].HV_Deflect_P,\
			para[i].HV_Cfilm,\
			para[i].HV_E,\
			para[i].HV_MCP1,\
			para[i].HV_MCP2,\
			para[i].Event_Num,\
			para[i].Event_Eff,\
			para[i].Event_Aux\
		);

		if (packet[i].ep2.stat){
			double q_dl = pow(2, -30);
			
			memcpy(&num32, &packet[i].ep2.flex[0], sizeof(uint8_t)*4);
			fill_Moon_J2000_q1 = int(B2L32(num32)) * q_dl; 
			para[i].Moon_J2000_q1 = fill_Moon_J2000_q1;

			memcpy(&num32, &packet[i].ep2.flex[4], sizeof(uint8_t)*4);
			fill_Moon_J2000_q2 = int(B2L32(num32)) * q_dl; 
			para[i].Moon_J2000_q2 = fill_Moon_J2000_q2;
			
			memcpy(&num32, &packet[i].ep2.flex[8], sizeof(uint8_t)*4);
			fill_Moon_J2000_q3 = int(B2L32(num32)) * q_dl; 
			para[i].Moon_J2000_q3 = fill_Moon_J2000_q3;

			memcpy(&num32, &packet[i].ep2.flex[12], sizeof(uint8_t)*4);
			fill_Moon_J2000_Star_X = int(B2L32(num32)); 
			para[i].Moon_J2000_Star_X = fill_Moon_J2000_Star_X;
			
			memcpy(&num32, &packet[i].ep2.flex[16], sizeof(uint8_t)*4);
			fill_Moon_J2000_Star_Y = int(B2L32(num32)); 
			para[i].Moon_J2000_Star_Y = fill_Moon_J2000_Star_Y;

			memcpy(&num32, &packet[i].ep2.flex[20], sizeof(uint8_t)*4);
			fill_Moon_J2000_Star_Z = int(B2L32(num32)); 
			para[i].Moon_J2000_Star_Z = fill_Moon_J2000_Star_Z;
			
			para[i].Star_Z_Solar_angle = fill_Star_Z_Solar_angle;
		}
		else{
			double unit_deg = 0.00001;
			memcpy(&num32, &packet[i].ep2.flex[18], sizeof(uint8_t)*4);
			fill_Star_Z_Solar_angle = int(B2L32(num32)) * unit_deg; 
			para[i].Star_Z_Solar_angle = fill_Star_Z_Solar_angle;
			
			para[i].Moon_J2000_q1 = fill_Moon_J2000_q1;
			para[i].Moon_J2000_q2 = fill_Moon_J2000_q2;
			para[i].Moon_J2000_q3 = fill_Moon_J2000_q3;
			para[i].Moon_J2000_Star_X = fill_Moon_J2000_Star_X;
			para[i].Moon_J2000_Star_Y = fill_Moon_J2000_Star_Y;
			para[i].Moon_J2000_Star_Z = fill_Moon_J2000_Star_Z;
		}

		//0x35
		if(packet[i].ep2.args[0] == 53){
			memcpy(&num32, &packet[i].ep2.args[1], sizeof(uint8_t)*4);
			fill_MCP1 = (int(B2L32(num32) & 0b111111111111)); 
			para[i].MCP1 = fill_MCP1;
		}
		else{
			para[i].MCP1 = fill_MCP1;

		}
		//0x37
		if(packet[i].ep2.args[0] == 55){
			memcpy(&num32, &packet[i].ep2.args[1], sizeof(uint8_t)*4);
			fill_AX = (int(B2L32(num32) & 0b111111111111)); 
			para[i].AX = fill_AX;
		}
		else{
			para[i].AX = fill_AX;
		}
		//0x39
		if(packet[i].ep2.args[0] == 57){
			memcpy(&num32, &packet[i].ep2.args[1], sizeof(uint8_t)*4);
			fill_MCP2 = (int(B2L32(num32) & 0b111111111111)); 
			para[i].MCP2 = fill_MCP2;
		}
		else{
			para[i].MCP2 = fill_MCP2;
		}
		//0x3B
	    	if(packet[i].ep2.args[0] == 59){
			memcpy(&num32, &packet[i].ep2.args[1], sizeof(uint8_t)*4);
			fill_AY = (int(B2L32(num32) & 0b111111111111)); 
			para[i].AY = fill_AY;
		}
		else{
			para[i].AY = fill_AY;
		}
		if(packet[i].ep2.args[0] == 142){
			fill_BX1 = signed_7bits(packet[i].ep2.args[2]); 
			para[i].BX1 = fill_BX1;
			fill_BX2 = signed_7bits(packet[i].ep2.args[4]); 
			para[i].BX2 = fill_BY1;
		}
		else{
			para[i].BX1 = fill_BX1;
			para[i].BX2 = fill_BX2;
		}
		if(packet[i].ep2.args[0] == 146){
			fill_BY1 = signed_7bits(packet[i].ep2.args[2]); 
			para[i].BY1 = fill_BY1;
			fill_BY2 = signed_7bits(packet[i].ep2.args[4]); 
			para[i].BY2 = fill_BY2;
		}
		else{
			para[i].BY1 = fill_BY1;
			para[i].BY2 = fill_BY2;
		}
		if(packet[i].ep2.args[0] == 150){
			fill_BG = int(packet[i].ep2.args[1]); 
			para[i].BG = fill_BG;
		}
		else
		{
			para[i].BG = fill_BG;
		}
	}
	//long len_id = std::min(idx_1.size(),idx_2.size(),idx_3.size(),idx_4.size());
	//auto len_id = std::min({idx_1.size(), idx_2.size(), idx_3.size(), idx_4.size()});
	//auto len_id2 = std::min({idx2_1.size(), idx2_2.size(), idx2_3.size(), idx2_4.size()});
	return {para, num_patch};
}



data_ary<GenaL0> raw_2_L0(Packet *packet, unsigned long long num_patch, unsigned long long max_packets)
{

	GenaL0 *dataL0 = (GenaL0 *)malloc(num_patch * sizeof(GenaL0));
	#pragma
	for (long i = 0; i < num_patch; i++)
	{
		parsePacket(packet[i], dataL0[i]);
	}

	// if last, use i-1's
	for (long i = 0; i < num_patch - 1; ++i)
	{
		// if true, this packet does not have antenna pointing info
		// use next packet's
		// if still 0 use next ...
		long next = i;
		// make sure utc time is initialized with 0
		while (dataL0[i].eps.utc_time == 0)
		{
			next = next + 1;
			dataL0[i].eps.utc_time = dataL0[next].eps.utc_time;
			dataL0[i].eps.ant_angle = dataL0[next].eps.ant_angle;
			dataL0[i].eps.ant_azimuth = dataL0[next].eps.ant_azimuth;
			dataL0[i].eps.Star_angular_momentum = dataL0[next].eps.Star_angular_momentum;
			dataL0[i].eps.Star_Z_Solar_angle = dataL0[next].eps.Star_Z_Solar_angle;
			dataL0[i].eps.Shadow = dataL0[next].eps.Shadow;
			if (next == max_packets)
			{
				break;
			}
		}
		next = i;
		while (dataL0[i].eps.Moon_J2000_q1 == 0)
		{
			next = next + 1;
			dataL0[i].eps.Moon_J2000_q1 = dataL0[next].eps.Moon_J2000_q1;
			dataL0[i].eps.Moon_J2000_q2 = dataL0[next].eps.Moon_J2000_q2;
			dataL0[i].eps.Moon_J2000_q3 = dataL0[next].eps.Moon_J2000_q3;
			dataL0[i].eps.Moon_J2000_Star_X = dataL0[next].eps.Moon_J2000_Star_X;
			dataL0[i].eps.Moon_J2000_Star_Y = dataL0[next].eps.Moon_J2000_Star_Y;
			dataL0[i].eps.Moon_J2000_Star_Z = dataL0[next].eps.Moon_J2000_Star_Z;
			if (next == max_packets)
			{
				break;
			}
		}
	}

	// add fill mark on L0 data
	for (int i = 0; i < num_patch; ++i)
	{
		for (int j = 0; j < NUM_EVENT_PACKETS; ++j)
		{
			// dump all satellite filled 0
			// 0 00b:connect
			// 1 01b:head
			// 2 10b:tail
			// 3 11b:indept
			bool indept_tail = dataL0[i].BX == 2 || dataL0[i].BX == 3;
			// connect and next one is tail;
			// if next is end of packet, treat it as true
			// what if tail packet is missing, next one would be head
			bool last_connect = dataL0[i].BX == 0 && (i + 1 < num_patch) ? (dataL0[i + 1].BX == 2) : (true);
			int DataSum = 0;
			for (int q = 0; q < 2; q++)
			{
				DataSum = DataSum + dataL0[i].evts[j].c1[q] + dataL0[i].evts[j].c1[q] + dataL0[i].evts[j].c5[q] + dataL0[i].evts[j].c7[q];
			}
			for (int q = 0; q < 4; q++)
			{
				DataSum = DataSum + dataL0[i].evts[j].tof[q];
			}
			DataSum = DataSum + dataL0[i].evts[j].ep;
			bool all_zero = DataSum == 0;
			if ((indept_tail || last_connect) && all_zero)
			{
				dataL0[i].evts[j].isfill = 1;
			}
			else
			{
				dataL0[i].evts[j].isfill = 0;
			}
		}
	}
	return {dataL0, num_patch};
}

data_ary<GenaL1> L0_2_L1(data_ary<GenaL0> ary)
{

	GenaL0 *l0 = ary.data;
	unsigned long num = ary.count;
	// nsigned long size_L1 = num * NUM_EVENT_PACKETS;
	unsigned long size_L1 = 0;
	vec_2int eff_ij;

	for (int i = 0; i < num; ++i)
	{
		for (int j = 0; j < NUM_EVENT_PACKETS; ++j)
		{
			if (l0[i].evts[j].isfill == 0)
			{
				eff_ij.push_back(make_pair(i, j));
				size_L1 = size_L1 + 1;
			}
		}
	}

	GenaL1 *l1 = (GenaL1 *)malloc(size_L1 * sizeof(GenaL1));
	for (int k = 0; k < size_L1; ++k)
	{
		int i = eff_ij[k].first;
		int j = eff_ij[k].second;
		l1[k].Tim_s = l0[i].Tim_s;
		l1[k].Tim_ms = l0[i].Tim_ms;
		l1[k].eps = l0[i].eps;
		l1[k].ep_mv = eplhe2ep(l0[i].evts[j].ep);
		for (int k = 0; k < 4; ++k)
		{
			l1[k].fly_time[k] = tdc2tof(l0[i].evts[j].tof[k]);
		}
        if(SN_ID == 1){
    		
if (l0[i].evts[j].prob_flag == 0)
		{
			l1[k].flag = 1;
			for (int q = 0; q < 2; ++q)
			{
				l1[k].tx1[q] = tdc2tof(l0[i].evts[j].c5[q]);
				l1[k].tx2[q] = tdc2tof(l0[i].evts[j].c7[q]);
				l1[k].ty1[q] = tdc2tof(l0[i].evts[j].c1[q]);
				l1[k].ty2[q] = tdc2tof(l0[i].evts[j].c3[q]);
			}
		}
		else if (l0[i].evts[j].prob_flag == 1)
		{
			l1[k].flag = 2;
			for (int q = 0; q < 2; ++q)
			{
				l1[k].tx1[q] = tdc2tof(l0[i].evts[j].c1[q]);
				l1[k].tx2[q] = tdc2tof(l0[i].evts[j].c3[q]);
				l1[k].ty1[q] = tdc2tof(l0[i].evts[j].c5[q]);
				l1[k].ty2[q] = tdc2tof(l0[i].evts[j].c7[q]);
			}
		}
}
else if (SN_ID == 2){
	
if (l0[i].evts[j].prob_flag == 0)
		{
			l1[k].flag = 1;
			for (int q = 0; q < 2; ++q)
			{
				l1[k].tx1[q] = tdc2tof(l0[i].evts[j].c5[q]);
				l1[k].tx2[q] = tdc2tof(l0[i].evts[j].c7[q]);
				l1[k].ty1[q] = tdc2tof(l0[i].evts[j].c1[q]);
				l1[k].ty2[q] = tdc2tof(l0[i].evts[j].c3[q]);
			}
		}
		else if (l0[i].evts[j].prob_flag == 1)
		{
			l1[k].flag = 2;
			for (int q = 0; q < 2; ++q)
			{
				l1[k].tx1[q] = tdc2tof(l0[i].evts[j].c7[q]);
				l1[k].tx2[q] = tdc2tof(l0[i].evts[j].c5[q]);
				l1[k].ty1[q] = tdc2tof(l0[i].evts[j].c3[q]);
				l1[k].ty2[q] = tdc2tof(l0[i].evts[j].c1[q]);
			}
		}

}
	}
	return {l1, size_L1};
}

bool rule_prob_a(const GenaL1 &l1)
{
	return l1.flag == 1;
}

bool rule_prob_b(const GenaL1 &l1)
{
	return l1.flag == 2;
}

// function on L1
// 1A data: anode 8 channel not all zero
bool rule_L1A(const GenaL1 &l1)
{
	double sum_anode = 0;
	for (int j = 0; j < 2; ++j)
	{
		sum_anode = sum_anode + l1.tx1[j] + l1.tx2[j] + l1.ty1[j] + l1.ty2[j];
	}
	return sum_anode != 0;
}

//1B data: all four start signal  > 0, and within limits
bool rule_L1B(const GenaL1 &l1)
{
	bool start_eff = (l1.tx1[0] > 0 && l1.tx2[0] > 0 && l1.ty1[0] > 0 && l1.ty2[0] > 0);
	//double ns_xlim;
	//double ns_ylim;
	//(l1.flag == 1) ? (ns_xlim = A_NS_XLIM) : (ns_xlim = B_NS_XLIM);
	//(l1.flag == 1) ? (ns_ylim = A_NS_YLIM) : (ns_ylim = B_NS_YLIM);
	//// start signal should less than limited
	//bool ns_eff = (l1.tx1[0] < ns_xlim && l1.tx2[0] < ns_xlim && l1.ty1[0] < ns_ylim && l1.ty2[0] < ns_ylim);
	//return start_eff && ns_eff;
	return start_eff;
}

//1C: based on 1B with all four stop nonzero
bool rule_L1C(const GenaL1 &l1)
{
	bool bool_1B = rule_L1B(l1);
	//bool start_eff = (l1.tx1[0] > 0 && l1.tx2[0] > 0 && l1.ty1[0] > 0 && l1.ty2[0] > 0);
	bool stop_eff = (l1.tx1[1] > 0 && l1.tx2[1] > 0 && l1.ty1[1] > 0 && l1.ty2[1] > 0);
	//return start_eff && stop_eff;
	return bool_1B && stop_eff;
}

// rational data used for generating L2 level data
bool rule_L1D(const GenaL1 &l1)
{

	// 1. based on 1C selection
	bool bool_1C = rule_L1C(l1);
	// 2. start time pair matched
	double x_sta = l1.tx1[0] + l1.tx2[0];
	double y_sta = l1.ty1[0] + l1.ty2[0];
	bool x_match = x_sta < 55 && x_sta > 45;
	bool y_match = y_sta < 45 && y_sta > 35;
	return bool_1C && x_match && y_match;
}

// 2A event: X1 Y1 X2 Y2 (mm) inside physical board
bool rule_L2A(const GenaL2 &l2)
{
    //bool first_eff = (abs(l2.x1) < X_CF / 2 && abs(l2.y1) < Y_CF / 2);
    bool first_eff = (abs(l2.x1) < X_MCP / 2 && abs(l2.y1) < Y_MCP / 2);
    bool second_eff = (abs(l2.x2) < X_MCP / 2 && abs(l2.y2) < Y_MCP / 2);
	
	//remove energy nan
	bool not_nan = !std::isnan(l2.energy);
	return first_eff && second_eff && not_nan;
}

// 2A event: X1 Y1 X2 Y2 (mm) inside physical board
bool rule_L2B(const GenaL2 &l2)
{
	bool bool_2A = rule_L2A(l2);
	bool energy_eff = l2.energy > 0 && l2.energy < 300 *1000;
	return bool_2A && energy_eff;
}


// L2 data
// convert to physical quantities
// time in ns
double ep_flytime(double V_CF, double V_MCP)
{
	// ev
	const double me = 511000;
	// m
	const double l_acc = 3.1 * 1e-3;
	const double l_fly = 47 * 1e-3;
	// 1
	double energy = abs(V_CF - V_MCP);
	// m/s
	double ve = sqrt(2 * energy / me) * SPEED_C;
	// ns
	double time_acc = 2 * l_acc / ve * 1e9;
	double time_fly = l_fly / ve * 1e9;
	return time_acc + time_fly;
}

void which_spec(GenaL2 &l2, double tof, double ep)
{
	int spec;
	double mH = 1836 * 511000.;
	double mO = 15.872162485065711 * mH;
	if (tof < 40)
	{
		spec = 0;
	}
	else if (tof > 66.3)
	{
		spec = 1;
	}
	else
	{
		if (ep < 1500)
		{
			spec = 0;
		}
		else if (ep > 1580)
		{
			spec = 1;
		}
		else
		{
			spec = 0;
		}
	}
	l2.spec = spec;
	(spec == 0) ? (l2.mass = mH) : (l2.mass = mO);
}

data_ary<GenaL2> L1_2_L2(data_ary<GenaL1> ary)
{
	GenaL1 *l1 = ary.data;
	unsigned long num = ary.count;
	GenaL2 *l2 = (GenaL2 *)malloc(num * sizeof(GenaL2));
	double V_Cf;
	double V_MCP;
	double Tof_delt_e;
	double X1;
	double Y1;
	double X2;
	double Y2;
	double tof_mcp1;
	double tof_anode_x;
	double tof_anode_y;
	double ena_tof;
	double ena_speed;
	double fly_dist, vx, vy;
	double sum0_tx1_tx2, sum0_ty1_ty2, sum1_tx1_tx2, sum1_ty1_ty2;
	for (long i = 0; i < num; ++i)
	{
		l2[i].flag = l1[i].flag;
		l2[i].Tim_s = l1[i].Tim_s;
		l2[i].Tim_ms = l1[i].Tim_ms;
		l2[i].Star_Z_Solar_angle = l1[i].eps.Star_Z_Solar_angle;
		l2[i].Moon_J2000_q1 = l1[i].eps.Moon_J2000_q1;
		l2[i].Moon_J2000_q2 = l1[i].eps.Moon_J2000_q2;
		l2[i].Moon_J2000_q3 = l1[i].eps.Moon_J2000_q3;

		l2[i].Moon_J2000_Star_X = l1[i].eps.Moon_J2000_Star_X;
		l2[i].Moon_J2000_Star_Y = l1[i].eps.Moon_J2000_Star_Y;
		l2[i].Moon_J2000_Star_Z = l1[i].eps.Moon_J2000_Star_Z;

		V_Cf = l1[i].eps.HV_Cfilm;
		V_MCP = l1[i].eps.HV_E;

		Tof_delt_e = ep_flytime(V_Cf, V_MCP);
		tof_mcp1 = l1[i].fly_time[0];

		X1 = l1[i].tx2[0] - l1[i].tx1[0];
		Y1 = l1[i].ty2[0] - l1[i].ty1[0];
		X2 = l1[i].tx2[1] - l1[i].tx1[1];
		Y2 = l1[i].ty2[1] - l1[i].ty1[1];

		sum0_tx1_tx2 = l1[i].tx1[0] + l1[i].tx2[0];
		sum0_ty1_ty2 = l1[i].ty1[0] + l1[i].ty2[0];
		sum1_tx1_tx2 = l1[i].tx1[1] + l1[i].tx2[1];
		sum1_ty1_ty2 = l1[i].ty1[1] + l1[i].ty2[1];

		// contains x sub and y sub, thus divided by 2
		tof_anode_x = (sum1_tx1_tx2 - sum0_tx1_tx2) / 2;
		tof_anode_y = (sum1_ty1_ty2 - sum0_ty1_ty2) / 2;
		if (l1[i].flag == 1)
		{
                    
			l2[i].x1 = (X1 - A_X_CENTER_NS) * A_X_COF;
			l2[i].y1 = (Y1 - A_Y_CENTER_NS) * A_Y_COF;
			l2[i].x2 = (X2 - A_X_CENTER_NS) * A_X_COF;
			l2[i].y2 = (Y2 - A_Y_CENTER_NS) * A_Y_COF;
		}
		if (l1[i].flag == 2)
		{
			l2[i].x1 = (X1 - B_X_CENTER_NS) * B_X_COF;
			l2[i].y1 = (Y1 - B_Y_CENTER_NS) * B_Y_COF;
			l2[i].x2 = (X2 - B_X_CENTER_NS) * B_X_COF;
			l2[i].y2 = (Y2 - B_Y_CENTER_NS) * B_Y_COF;
		}
		l2[i].ep_mv = l1[i].ep_mv;
		// discuss fly time;
		fly_dist = sqrt(pow((l2[i].x1 - l2[i].x2), 2) + pow((l2[i].y1 - l2[i].y2), 2) + pow(50, 2));
		//cout<<tof_mcp1<<endl;
		//most are zero
		if (tof_mcp1 > 0 && abs(tof_anode_x - tof_anode_y) <= 5 && l1[i].tx1[1] > tof_mcp1 && l1[i].tx2[1] > tof_mcp1 && l1[i].ty1[1] > tof_mcp1 && l1[i].ty2[1] > tof_mcp1)
		{
			ena_tof = tof_mcp1;
			which_spec(l2[i], ena_tof, l2[i].ep_mv);
			// if tof are close to tof obtained from annode data, use tof
			if (abs(tof_mcp1 - tof_anode_x) < 5 and abs(tof_mcp1 - tof_anode_y) < 5)
			{
				ena_tof = tof_mcp1 + Tof_delt_e;
				ena_speed = (fly_dist * 1e-3) / (ena_tof * 1e-9) / SPEED_C;
				l2[i].fly_time = ena_tof;
				l2[i].energy = (l2[i].mass * pow(ena_speed, 2) / 2);
			}
			// else tof is not match with annode, use annode (max one)
			else
			{
				ena_tof = (tof_anode_x > tof_anode_y) ? (tof_anode_y) : (tof_anode_x) + Tof_delt_e;
				which_spec(l2[i], ena_tof, l2[i].ep_mv);
				ena_speed = (fly_dist * 1e-3) / (ena_tof * 1e-9) / SPEED_C;
				l2[i].fly_time = ena_tof;
				l2[i].energy = (l2[i].mass * pow(ena_speed, 2) / 2);
			}
		}
		// if no tof or have tof but > 5
		else if (abs(tof_anode_x - tof_anode_y) < 8)
		{
			// use average anode tof
			ena_tof = (tof_anode_x + tof_anode_y) / 2;
			//which_spec(l2[i], tof_mcp1, l2[i].ep_mv);
			which_spec(l2[i], ena_tof, l2[i].ep_mv);
			ena_tof = (tof_anode_x > tof_anode_y) ? (tof_anode_y) : (tof_anode_x) + Tof_delt_e;
			ena_speed = (fly_dist * 1e-3) / (ena_tof * 1e-9) / SPEED_C;
			l2[i].fly_time = ena_tof;
			l2[i].energy = (l2[i].mass * pow(ena_speed, 2) / 2);
		}
		else
		{
			l2[i].mass = 0;
			l2[i].fly_time = NAN;
			l2[i].energy = NAN;
		}
			//cout << l2[i].energy << endl;
	}
	return {l2, num};
}
