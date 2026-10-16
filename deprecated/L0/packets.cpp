#include <cstdint>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "structures.hpp"
using namespace std;

uint16_t B2L16(uint16_t num) {
	return ((num >> 8) & 0xff) | ((num << 8) & 0xff00);
}

uint32_t B2L32(uint32_t num) {
	return ((num >> 24) & 0xff) |
		((num >> 8) & 0xff00) |
		((num << 8) & 0xff0000) |
		((num << 24) & 0xff000000);
}

uint64_t B2L64(uint64_t num) {
	return ((num >> 56) & 0xff) |
		((num >> 40) & 0xff00) |
		((num >> 24) & 0xff0000) |
		((num >> 8) & 0xff000000) |
		((num << 8) & 0xff00000000) |
		((num << 24) & 0xff0000000000) |
		((num << 40) & 0xff000000000000) |
		((num << 56) & 0xff00000000000000);
}

int signed_7bits(uint8_t bits){
	int sign = (bits & 0x40) ? -1 : 1;
	int value = bits & 0x3F;
	return sign * value;
}
void seek_ep1(Ep1 & ep, double & HV_Deflect_N, double & HV_Deflect_P, double & HV_Cfilm, double & HV_E, double & HV_MCP1, double & HV_MCP2,\
		unsigned short & Event_Num, unsigned short & Event_Eff, unsigned short & Event_Aux) 
{
	double const1 = 50000./32768;
	double const2 = 10000./32768;
	// 偏转负高压：补码，[-15000,0]
	HV_Deflect_N = short(B2L16((ep.HVs[0]))) * const1;
	// 偏转正高压：[0,15000]
	HV_Deflect_P = (unsigned short)(B2L16((ep.HVs[1]))) *  const1;
	// 碳膜主高压：补码，[-5000,0]
	HV_Cfilm = short(B2L16((ep.HVs[2]))) * const2;
	// 引导电场高压：补码，[-3000,0]
	HV_E = short(B2L16((ep.HVs[3]))) * const2;
	// MCP1高压：补码，[-3000,0]
	HV_MCP1 = short(B2L16((ep.HVs[4]))) * const2;
	// MCP2高压：补码，[-3000,0]
	HV_MCP2 =  short(B2L16((ep.HVs[5]))) * const2;
	// 假设events_counts中的每个元素是2字节（根据实际调整类型和长度
	Event_Num = B2L16(ep.event_count[0]);
	Event_Eff = B2L16(ep.event_count[1]);
	Event_Aux = B2L16(ep.event_count[2]);
}

void seek_ep2(Ep2 &ep, double & utc_time, double & ant_angle, double & ant_azimuth, double & Star_angular_momentum ,double &  Star_Z_Solar_angle ,unsigned int & Shadow ,\
		double & Moon_J2000_q1 ,double & Moon_J2000_q2 ,double & Moon_J2000_q3 ,int & Moon_J2000_Star_X ,int & Moon_J2000_Star_Y ,int & Moon_J2000_Star_Z)
{
	if (ep.stat){
		double q_dl = pow(2, -30);
		uint32_t num32;

		memcpy(&num32, &ep.flex[0], sizeof(uint8_t)*4);
		Moon_J2000_q1 = int(B2L32(num32)) * q_dl; 

		memcpy(&num32, &ep.flex[4], sizeof(uint8_t)*4);
		Moon_J2000_q2 = int(B2L32(num32))  * q_dl; 

		memcpy(&num32, &ep.flex[8], sizeof(uint8_t)*4);
		Moon_J2000_q3 = int(B2L32(num32))  * q_dl; 

		memcpy(&num32, &ep.flex[12], sizeof(uint8_t)*4);
		Moon_J2000_Star_X = int(B2L32(num32)); 

		memcpy(&num32, &ep.flex[16], sizeof(uint8_t)*4);
		Moon_J2000_Star_Y = int(B2L32(num32)); 

		memcpy(&num32, &ep.flex[20], sizeof(uint8_t)*4);
		Moon_J2000_Star_Z =  int(B2L32(num32)); 
	
		utc_time = 0;
		ant_angle = 0;
		ant_azimuth = 0;
		Star_angular_momentum = 0; 
		Star_Z_Solar_angle = 0;
		Shadow = 0;
	}
	else{
		uint64_t num64;
		memcpy(&num64, &ep.flex[0], sizeof(uint8_t)*6);
		utc_time = B2L64(num64 << (8 * (8-6))) *0.01;

		double unit_deg = 0.00001;
		uint32_t num32;
		memcpy(&num32, &ep.flex[6], sizeof(uint8_t)*4);
		ant_angle = int(B2L32(num32)) * unit_deg; 

		memcpy(&num32, &ep.flex[10], sizeof(uint8_t)*4);
		ant_azimuth = int(B2L32(num32)) * unit_deg; 

		memcpy(&num32, &ep.flex[14], sizeof(uint8_t)*4);
		Star_angular_momentum = int(B2L32(num32)) * unit_deg; 

		memcpy(&num32, &ep.flex[18], sizeof(uint8_t)*4);
		Star_Z_Solar_angle = int(B2L32(num32)) * unit_deg; 

		Shadow = ep.flex[22];
		Moon_J2000_q1 = 0; 
		Moon_J2000_q2 = 0; 
		Moon_J2000_q3 = 0; 
		Moon_J2000_Star_X = 0; 
		Moon_J2000_Star_Y = 0; 
		Moon_J2000_Star_Z = 0; 
	}
}

void seek_event(Event& ev, unsigned int& event_D_mark, unsigned int& event_Incid,unsigned int& event_D_incid_m, unsigned int& event_D_incid_z_m, unsigned  int& event_Times, unsigned  int& event_Timeus, unsigned int& event_EplHe, unsigned int event_fly_time[][2],unsigned  int event01_data[][2],unsigned  int event03_data[][2],unsigned  int event05_data[][2], unsigned int event07_data[][2]){
	// prob flag，0:A 1:B 取标志位，高位第1位
	event_D_mark = ev.flag >> 7;
	// 事件有效性，取标志位，高位第2位
	event_Incid = (ev.flag >> 6) & 1;
	// 探测器探测模式标记，取标志位，高位第3、4位
	event_D_incid_m = (ev.flag >> 4) & 0b11;
	// 探测器探测子模式标记，低4位
	event_D_incid_z_m = ev.flag & 0b1111;

	//internal time code; 
	event_Times = int(B2L16(ev.time_ms));
	event_Timeus = int(B2L16(ev.time_us));
	// connect three bytes
	event_EplHe = ((uint32_t)ev.se_amp[0] << 16) | ((uint32_t)ev.se_amp[1] << 8) | (uint32_t)ev.se_amp[2];
	// 进行位与操作，提取相应位
	event_EplHe &= 0b111111111111;  

	//use array[*][2] to store time
	//[*][0] represents group number
	//[*][1] represents tdc time

	// tof0 fly time measurements
	// 4 groups identical，
	unsigned int ch_fly_time[4][2];

	// mcp data, 
	//1 3 5 7: represnets 4 sides
	// each side have 2 group: 0 denotes second electron tdc, 1 denotes ENA tdc
	unsigned int ch01_data[2][2];
	unsigned int ch03_data[2][2];
	unsigned int ch05_data[2][2];
	unsigned int ch07_data[2][2];

	//channel 0
	uint16_t ev_c;
	//ii is channel index
	for (int i = 0; i < 4; ++i) {
		ev_c = B2L16(ev.c0[i]);
		ch_fly_time[i][0] = ev_c >> 13;
		ch_fly_time[i][1] = ev_c & 0b1111111111111;
	}
	//channel 1
	//ij ... is group index
	for (int i = 0; i < 2; ++i) {
		ev_c = B2L16(ev.c1[i]);
		ch01_data[i][0] = ev_c >> 13;
		ch01_data[i][1] = ev_c & 0b1111111111111;
	}
	//channel 3
	for (int i = 0; i < 2; ++i) {
		ev_c = B2L16(ev.c3[i]);
		ch03_data[i][0] = ev_c >> 13;
		ch03_data[i][1] = ev_c & 0b1111111111111;
	}
	//channel 5
	for (int i = 0; i < 2; ++i) {
		ev_c =  B2L16(ev.c5[i]);
		ch05_data[i][0] = ev_c >> 13;
		ch05_data[i][1] = ev_c & 0b1111111111111;
	}
	//channel 7
	for (int i = 0; i < 2; ++i) {
		ev_c = B2L16(ev.c7[i]);
		ch07_data[i][0] = ev_c >> 13;
		ch07_data[i][1] = ev_c & 0b1111111111111;
	}

	memcpy(event_fly_time, ch_fly_time, sizeof(ch_fly_time));
	memcpy(event01_data, ch01_data, sizeof(ch01_data));
	memcpy(event03_data, ch03_data, sizeof(ch03_data));
	memcpy(event05_data, ch05_data, sizeof(ch05_data));
	memcpy(event07_data, ch07_data, sizeof(ch07_data));
}



//In packet 
//Out L0
void parsePacket(Packet& packet, GenaL0& l0) {
	uint16_t seq0;
	//group mark
	//01b first packet
	//00b connect packet
	//10b last packet
	//11b indepacket packet?
	seq0 = B2L16(packet.seq[0]);
	l0.BX = seq0 >> 14 & 0b11;
	//packet count
	l0.BM = seq0 & 0b11111111111111;
	l0.Tim_s = B2L32(packet.time_s);
	l0.Tim_ms = B2L16(packet.time_ms);
	seek_ep1(packet.ep1,\
			l0.eps.HV_Deflect_N,\
			l0.eps.HV_Deflect_P,\
			l0.eps.HV_Cfilm,\
			l0.eps.HV_E,\
			l0.eps.HV_MCP1,\
			l0.eps.HV_MCP2,\
			l0.eps.Event_Num,\
			l0.eps.Event_Eff,\
			l0.eps.Event_Aux\
		);
	seek_ep2(packet.ep2,\
			l0.eps.utc_time,\
			l0.eps.ant_angle,\
			l0.eps.ant_azimuth,\
			l0.eps.Star_angular_momentum,\
			l0.eps.Star_Z_Solar_angle,\
			l0.eps.Shadow,\
			l0.eps.Moon_J2000_q1,\
			l0.eps.Moon_J2000_q2,\
			l0.eps.Moon_J2000_q3,\
			l0.eps.Moon_J2000_Star_X,\
			l0.eps.Moon_J2000_Star_Y,\
			l0.eps.Moon_J2000_Star_Z\
		);

	unsigned int event_D_mark,event_Incid,event_D_incid_m,event_D_incid_z_m;
	unsigned int event_Timems,event_Timeus;
	unsigned int event_EplHe;
	unsigned int event_fly_time[4][2];
	// 1 3 5 7 each have 2 group and 2 channel
	//dim1: channel number, tdc time
	//dim2: start, end
	unsigned int event01_data[2][2];
	unsigned int event03_data[2][2];
	unsigned int event05_data[2][2];
	unsigned int event07_data[2][2];
	for(int i=0;i<NUM_EVENT_PACKETS;++i){
		//seek_event()
		seek_event(packet.evts[i],event_D_mark,event_Incid,event_D_incid_m,event_D_incid_z_m,\
				event_Timems,event_Timeus,event_EplHe, event_fly_time, event01_data, event03_data, event05_data, event07_data);
		l0.evts[i].prob_flag = event_D_mark;
		l0.evts[i].eff_flag = event_Incid;
		l0.evts[i].ABmod= event_D_incid_m;
		l0.evts[i].submod= event_D_incid_z_m;
		l0.evts[i].ep = event_EplHe;
		l0.evts[i].internal_time = event_Timems + event_Timeus*1e-6;
		for(int j=0; j<4; ++j){
			l0.evts[i].tof[j] = event_fly_time[j][1];
		}
		for(int j=0;j<2;++j){
			l0.evts[i].c1[j] = event01_data[j][1];
			l0.evts[i].c3[j] = event03_data[j][1];
			l0.evts[i].c5[j] = event05_data[j][1];
			l0.evts[i].c7[j] = event07_data[j][1];
		}
	}
};


//vector<int64_t> get_packet_header_addr(unsigned int * buffer, long len_range[2]) {
//Already in memory
vector<int64_t> get_packet_header_addr(uint8_t * buffer, long long len_range[2]) {
	vector<int64_t> headseek_num;
	//fseek(f,len_range[0], SEEK_SET);
	long remaining_bytes = len_range[1] -  len_range[0];
	long cur=0;
	unsigned int headdata[4];
	//int count=0;
	long read_size;
	while (remaining_bytes > 0) {
		//count++;
		read_size = remaining_bytes < 4? remaining_bytes : 4;
		for(int i=0;i<read_size;i++){
			//memcpy(headdata, buffer + cur, read_size*sizeof(uint8_t));
			headdata[i] = buffer[cur+i];
		}
		cur = cur + read_size;
		if (read_size < 4) {
			break;
		}
		// check if all 4 bytes readed match E2 25 0C 4E
		if (headdata[0] == 226 && headdata[1] == 37 && headdata[2] == 12 && headdata[3] == 78) {
			headseek_num.push_back(cur - 4);
			continue;
		}
		// if missing three bytes
		else if (headdata[1] == 226 && headdata[2] == 37 && headdata[3] == 12) {
			cur = cur - 3;
		}
		// same, if missing two bytes
		else if (headdata[2] == 226 && headdata[3] == 37) {
			cur = cur -2;
		}   
		// however, if there is missing one Byte, then the last one of the 4 B would be E5, which belong to next 2064
		else if (headdata[3] == 226) {
			cur = cur -1;
		}
		remaining_bytes = len_range[1] - len_range[0] - cur;
	}
	//cout<< count <<endl;
	return headseek_num;
}

vector<int64_t> find_2064(const vector<int64_t>& a, int64_t L) {
	vector<int64_t> result;
	if (a.size() > 0) {
		for (long i = 1; i < a.size(); ++i) {
			if (a[i] - a[i - 1] == 2064) {
				result.push_back(a[i - 1]);
			}
		}
		// 判断最后一个元素是否满足条件（文件总长度L与最后一个元素的间隔是否为2064）
		if (L - a.back() == 2064) {
			result.push_back(a.back());
		}
	}
	return result;
}
