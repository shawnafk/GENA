#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include "structures.hpp"
#include <csv2/writer.hpp>
#include "io.hpp"

#include <iomanip>
#include <chrono>
#include <ctime>
using namespace csv2;
using namespace std;

std::string stamp2string(std::time_t seconds, int milliseconds)
{
	//1970-01-01 00:00:00 is at UTC timezone
	
	//Our readed seconds is to
	//TO BJ TIME 2021-01-01 00:00:00
	//TO UTC TIME IS 2020-12-31 16:00:00
	// 2021年1月1日00:00:00对应的BJ时间戳（以秒为单位）or 20201231160000 UTC TIMESTAMP
	//ANY WAY both have 1609430400
	//UTC 2021-1-1 00:00:00, BJ 2021-1-1 08:00:00; so more 8*3600 is added
	//std::time_t baseTime = 1609459200;
	std::time_t baseTime = 1609430400;
	// 计算出对应时间点的总秒数（加上传入的秒数和毫秒换算后的秒数）
	std::time_t targetSeconds = baseTime + seconds + milliseconds / 1000;
	std::chrono::system_clock::time_point tp = std::chrono::system_clock::from_time_t(targetSeconds);
	std::time_t updatedSeconds = std::chrono::system_clock::to_time_t(tp);
	std::tm *timeinfo = std::gmtime(&updatedSeconds);
	if (timeinfo == nullptr)
	{
		return "";
	}
	std::stringstream ss;
	// 获取年、月、日、时、分、秒的值
	int year = timeinfo->tm_year + 1900;
	int month = timeinfo->tm_mon + 1;
	int day = timeinfo->tm_mday;
	int hour = timeinfo->tm_hour;
	int minute = timeinfo->tm_min;
	int second = timeinfo->tm_sec;

	// 按照"YYYY-MM-DD HH:MM:SS.ms"格式构建字符串
	ss << std::setw(4) << std::setfill('0') << year << "-"
		<< std::setw(2) << std::setfill('0') << month << "-"
		<< std::setw(2) << std::setfill('0') << day << "-"
		<< std::setw(2) << std::setfill('0') << hour << ":"
		<< std::setw(2) << std::setfill('0') << minute << ":"
		<< std::setw(2) << std::setfill('0') << second << "."
		<< std::setw(3) << std::setfill('0') << milliseconds;

	return ss.str();
}

std::string ss_int(int num)
{
	std::stringstream ss;
	ss << num;
	std::string str_int = ss.str();
	return str_int;
}

std::string ss_double(double num)
{
	std::stringstream ss;
	ss.precision(6);
	ss << std::fixed << num;
	std::string str_double = ss.str();
	return str_double;
}

int Init_csv(const std::string &fileName, int level)
{
	std::vector<std::string> csv_header;
	std::ofstream stream(fileName);
	Writer<delimiter<','>> writer(stream);
	if (level == 0)
	{
		csv_header =
		{
			//"分组标识",
			"Group mark", 
			//"包序",
			"Packet order",          
			//"北京时间",
			"Time",          
			//"星时(s)",
			"Satellite seconds",     
			//"星时(ms)",
			"Satellite milseconds",  
			//"中性原子偏转负高压(V)",
			"HV_Deflect_N",          
			//"中性原子偏转正高压(V)",
			"HV_Deflect_P",          
			//"中性原子碳膜主高压(V)",
			"HV_Cfilm",              
			//"中性原子引导电场高压(V)",
			"HV_E",                  
			//"MCP1工作高压(V)",
			"HV_MCP1",               
			//"MCP2工作高压(V)",
			"HV_MCP1",               
			//"粒子事件数",
			"Event num",             
			//"有效事件数",
			"Event eff",             
			//"辅助事件数",
			"Event aux",             
			//"天线指向中心角(°)",
			"Ant_angle",             
			//"天线指向方位角(°)",
			"Ant_azimuth",           
			//"整星角动量计算值(Nms)",
			"Star_angular_momentum", 
			//"整星-Z轴与太阳矢量夹角(°)",
			"Star_Z_Solar_angle",    
			//"阴影判别标识(0/1)",
			"Shadow",                
			//"卫星本体机械坐标系相对于月心J2000惯性坐标系下的四元数q1",
			"Moon_J2000_q1",         
			//"卫星本体机械坐标系相对于月心J2000惯性坐标系下的四元数q2",
			"Moon_J2000_q2",         
			//"卫星本体机械坐标系相对于月心J2000惯性坐标系下的四元数q3",
			"Moon_J2000_q3",         
			//"月心J2000惯性系卫星X方向位置(m)",
			"Moon_J2000_Star_X",     
			//"月心J2000惯性系卫星Y方向位置(m)",
			"Moon_J2000_Star_Y",     
			//"月心J2000惯性系卫星Z方向位置(m)",
			"Moon_J2000_Star_Z",     
			//"探头标志位",
			"prob_flag",             
			//"有效辅助事件(1/0)",
			"eff_flag",              
			//"探头工作模式(A/B工作)",
			"ABmod",                 
			//"探头工作子模式",
			"Submod",                
			//"内部时间",
			"Internal time",
			"Tof_0",
			"Tof_1",
			"Tof_2",
			"Tof_3",
			"Ch01_0",
			"Ch01_1",
			"Ch03_0",
			"Ch03_1",
			"Ch05_0",
			"Ch05_1",
			"Ch07_0",
			"Ch07_1",
			"Secondary electron Amp",
			"IsFill"
			};
	}
	else if (level == 1)
	{
		csv_header = {
			"Date",
			"Timestamp",
			"HV_Deflect_N",
			"HV_Deflect_P",
			"HV_Cfilm",
			"HV_E",
			"HV_MCP1",
			"Ant_angle",
			"Ant_azimuth",
			"Star_angular_momentum",
			"Star_Z_Solar_angle",
			"Shadow",
			"Moon_J2000_q1",
			"Moon_J2000_q2",
			"Moon_J2000_q3",
			"Moon_J2000_Star_X",
			"Moon_J2000_Star_Y",
			"Moon_J2000_Star_Z",
			"Tof_0 (ns)",
			"Tof_1 (ns)",
			"Tof_2 (ns)",
			"Tof_3 (ns)",
			"X1_start (ns)",
			"X1_end (ns)",
			"X2_start (ns)",
			"X2_end (ns)",
			"Y1_start (ns)",
			"Y1_end (ns)",
			"Y2_start (ns)",
			"Y2_end (ns)",
			"Second electron amp (mv)",
			"ProbAB"
			};
	}
	else if (level == 2)
	{
		csv_header = {
			"Date",
			"Timestamp",
			"Star_Z_Solar_angle",
			"Moon_J2000_q1",
			"Moon_J2000_q2",
			"Moon_J2000_q3",
			"Moon_J2000_Star_X",
			"Moon_J2000_Star_Y",
			"Moon_J2000_Star_Z",
			"TOF (ns)",
			"X1 (mm)",
			"Y1 (mm)",
			"X2 (mm)",
			"Y2 (mm)",
			"Energy",
			"SPEC",
			"SE (mv)"};
	}
	writer.write_row(csv_header);
	stream.close();
	return 0;
}

int Init_csv_p(const std::string &fileName, int level)
{
	std::vector<std::string> csv_header;
	std::ofstream stream(fileName);
	Writer<delimiter<','>> writer(stream);
	if (level == 0)
	{
		csv_header =
		{
			"Date",
			"Timestamp",
			"Internal time",
			"Tof_0",
			"Tof_1",
			"Tof_2",
			"Tof_3",
			"Ch01_0",
			"Ch01_1",
			"Ch03_0",
			"Ch03_1",
			"Ch05_0",
			"Ch05_1",
			"Ch07_0",
			"Ch07_1",
			"EP",
			"IsFill"
		};
	}
	else if (level == 1)
	{
		csv_header = {
			"Date",
			"Timestamp",
			"Tof_0",
			"Tof_1",
			"Tof_2",
			"Tof_3",
		    "TX_1s",
			"TX_2s",
			"TY_1s",
			"TY_2s",
			"TX_1e",
			"TX_2e",
			"TY_1e",
			"TY_2e",
			"EP"
		};
	}
	else if (level == 2)
	{
		csv_header = {
			"Date",
			"Timestamp",
			"ToF",
			"X1",
			"Y1",
			"X2",
			"Y2",
			"Energy",
			"SPEC",
			"SE"
		};
	}
	writer.write_row(csv_header);
	stream.close();
	return 0;
}

int Init_csv_parameter(const std::string &fileName)
{
	std::vector<std::string> csv_header;
	std::ofstream stream(fileName);
	Writer<delimiter<','>> writer(stream);
	csv_header =
	{
			"Date",
			"Timestamp",
			"MCP1",
			"MCP2",
			"AX",
			"AY",
			"BX1",
			"BX2",
			"BY1",
			"BY2",
			"BG",
			"TA",
			"TB",
			"HV_Deflect_N",
			"HV_Deflect_P",
			"HV_Cfilm",
			"HV_E",
			"HV_MCP1",
			"HV_MCP2",
			"Event_Num",
			"Event_Eff",
			"Event_Aux",
			"Star_Z_Solar_angle",
			"Moon_J2000_q1",
			"Moon_J2000_q2",
			"Moon_J2000_q3",
			"Moon_J2000_Star_X",
			"Moon_J2000_Star_Y",
			"Moon_J2000_Star_Z",
	};
	writer.write_row(csv_header);
	stream.close();
	return 0;
}

int Output_parameter_csv(data_ary<GenaParameters> ary, const std::string &fileName)
{
	size_t num_events = ary.count;
	GenaParameters *para = ary.data;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);

	int columns = 29;
	std::vector<std::vector<std::string>> rows(num_events, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_events; i++)
	{
		rows[i][0] = stamp2string(int(para[i].Timestamp), 0);
		rows[i][1] = ss_double(para[i].Timestamp + 1609459200);
		rows[i][2] = ss_double(para[i].MCP1);
		rows[i][3] = ss_double(para[i].MCP2);
		rows[i][4] = ss_double(para[i].AX);
		rows[i][5] = ss_double(para[i].AY);
		rows[i][6] = ss_int(para[i].BX1);
		rows[i][7] = ss_int(para[i].BX2);
		rows[i][8] = ss_int(para[i].BY1);
		rows[i][9] = ss_int(para[i].BY2);
		rows[i][10] = ss_int(para[i].BG);
		rows[i][11] =	ss_double(para[i].TA);
		rows[i][12] =	ss_double(para[i].TB);
		rows[i][13] =	ss_double(para[i].HV_Deflect_N);
		rows[i][14] =	ss_double(para[i].HV_Deflect_P);
		rows[i][15] =	ss_double(para[i].HV_Cfilm);
		rows[i][16] =	ss_double(para[i].HV_E);
		rows[i][17] =	ss_double(para[i].HV_MCP1);
		rows[i][18] =	ss_double(para[i].HV_MCP2);
		rows[i][19] =	ss_double(para[i].Event_Num);
		rows[i][20] =	ss_double(para[i].Event_Eff);
		rows[i][21] =	ss_double(para[i].Event_Aux);
		rows[i][22] = ss_double(para[i].Star_Z_Solar_angle);
		rows[i][23] = ss_double(para[i].Moon_J2000_q1);
		rows[i][24] = ss_double(para[i].Moon_J2000_q2);
		rows[i][25] = ss_double(para[i].Moon_J2000_q3);
		rows[i][26] = ss_int(para[i].Moon_J2000_Star_X);
		rows[i][27] = ss_int(para[i].Moon_J2000_Star_Y);
		rows[i][28] = ss_int(para[i].Moon_J2000_Star_Z);
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}

int OutputL0_csv(data_ary<GenaL0> ary, const std::string &fileName)
{
	GenaL0 *l0 = ary.data;
	size_t num_packets = ary.count;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);
	int columns = 44;
	// 初始化二维vector，外层vector大小为num * 62，代表总共有num * 62行数据
	std::vector<std::vector<std::string>> rows(num_packets * 62, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_packets * 62; i++)
	{
		// 计算当前属于第几组数据（每组对应62行）
		size_t packet_idx = i / 62;
		// 计算在当前组内是第几行
		int row_in_packet = i % 62;
		rows[i][0] = ss_int(l0[packet_idx].BX);
		rows[i][1] = ss_int(l0[packet_idx].BM);
		// BJ time
		rows[i][2] = stamp2string(l0[packet_idx].Tim_s, l0[packet_idx].Tim_ms);
		rows[i][3] = ss_int(l0[packet_idx].Tim_s);
		rows[i][4] = ss_int(l0[packet_idx].Tim_ms);
		rows[i][5] = ss_double(l0[packet_idx].eps.HV_Deflect_N);
		rows[i][6] = ss_double(l0[packet_idx].eps.HV_Deflect_P);
		rows[i][7] = ss_double(l0[packet_idx].eps.HV_Cfilm);
		rows[i][8] = ss_double(l0[packet_idx].eps.HV_E);
		rows[i][9] = ss_double(l0[packet_idx].eps.HV_MCP1);
		rows[i][10] = ss_double(l0[packet_idx].eps.HV_MCP2);
		rows[i][11] = ss_int(l0[packet_idx].eps.Event_Num);
		rows[i][12] = ss_int(l0[packet_idx].eps.Event_Eff);
		rows[i][13] = ss_int(l0[packet_idx].eps.Event_Aux);
		rows[i][14] = ss_double(l0[packet_idx].eps.ant_angle);
		rows[i][15] = ss_double(l0[packet_idx].eps.ant_azimuth);
		rows[i][16] = ss_double(l0[packet_idx].eps.Star_angular_momentum);
		rows[i][17] = ss_double(l0[packet_idx].eps.Star_Z_Solar_angle);
		rows[i][18] = ss_int(l0[packet_idx].eps.Shadow);
		rows[i][19] = ss_double(l0[packet_idx].eps.Moon_J2000_q1);
		rows[i][20] = ss_double(l0[packet_idx].eps.Moon_J2000_q2);
		rows[i][21] = ss_double(l0[packet_idx].eps.Moon_J2000_q3);
		rows[i][22] = ss_int(l0[packet_idx].eps.Moon_J2000_Star_X);
		rows[i][23] = ss_int(l0[packet_idx].eps.Moon_J2000_Star_Y);
		rows[i][24] = ss_int(l0[packet_idx].eps.Moon_J2000_Star_Z);

		rows[i][25] = ss_int(l0[packet_idx].evts[row_in_packet].prob_flag);
		rows[i][26] = ss_int(l0[packet_idx].evts[row_in_packet].eff_flag);
		rows[i][27] = ss_int(l0[packet_idx].evts[row_in_packet].ABmod);
		rows[i][28] = ss_int(l0[packet_idx].evts[row_in_packet].submod);
		rows[i][29] = ss_double(l0[packet_idx].evts[row_in_packet].internal_time);

		for (int j = 0; j < 4; j++)
		{
			rows[i][30 + j] = ss_int(l0[packet_idx].evts[row_in_packet].tof[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			rows[i][34 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c1[j]);
			rows[i][36 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c3[j]);
			rows[i][38 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c5[j]);
			rows[i][40 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c7[j]);
		}
		rows[i][42] = ss_int(l0[packet_idx].evts[row_in_packet].ep);
		rows[i][43] = ss_int(l0[packet_idx].evts[row_in_packet].isfill);
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}

int OutputL1_csv(data_ary<GenaL1> ary, const std::string &fileName)
{
	GenaL1 *l1 = ary.data;
	size_t num_events = ary.count;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);

	int columns = 32;
	std::vector<std::vector<std::string>> rows(num_events, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_events; i++)
	{
		rows[i][0] = stamp2string(l1[i].Tim_s, l1[i].Tim_ms);
		rows[i][1] = ss_double(l1[i].Tim_s+l1[i].Tim_ms/1000. + 1609459200);
		rows[i][2] = ss_double(l1[i].eps.HV_Deflect_N);
		rows[i][3] = ss_double(l1[i].eps.HV_Deflect_P);
		rows[i][4] = ss_double(l1[i].eps.HV_Cfilm);
		rows[i][5] = ss_double(l1[i].eps.HV_E);
		rows[i][6] = ss_double(l1[i].eps.HV_MCP1);
		rows[i][7] = ss_double(l1[i].eps.ant_angle);
		rows[i][8] = ss_double(l1[i].eps.ant_azimuth);
		rows[i][9] = ss_double(l1[i].eps.Star_angular_momentum);
		rows[i][10] = ss_double(l1[i].eps.Star_Z_Solar_angle);
		rows[i][11] = ss_int(l1[i].eps.Shadow);
		rows[i][12] = ss_double(l1[i].eps.Moon_J2000_q1);
		rows[i][13] = ss_double(l1[i].eps.Moon_J2000_q2);
		rows[i][14] = ss_double(l1[i].eps.Moon_J2000_q3);
		rows[i][15] = ss_int(l1[i].eps.Moon_J2000_Star_X);
		rows[i][16] = ss_int(l1[i].eps.Moon_J2000_Star_Y);
		rows[i][17] = ss_int(l1[i].eps.Moon_J2000_Star_Z);

		for (int j = 0; j < 4; ++j)
		{
			rows[i][18 + j] = ss_double(l1[i].fly_time[j]);
		}
		for (int j = 0; j < 2; ++j)
		{
			// std::cout<<l1[i].tx1[j]<<std::endl;
			rows[i][22 + j] = ss_double(l1[i].tx1[j]);
			rows[i][24 + j] = ss_double(l1[i].tx2[j]);
			rows[i][26 + j] = ss_double(l1[i].ty1[j]);
			rows[i][28 + j] = ss_double(l1[i].ty2[j]);
		}
		rows[i][30] = ss_double(l1[i].ep_mv);
		rows[i][31] = (l1[i].flag==1)?("A"):("B");
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}

int OutputL2_csv(data_ary<GenaL2> ary, const std::string &fileName)
{
	GenaL2 *l2 = ary.data;
	size_t num_events = ary.count;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);

	int columns = 17;
	std::vector<std::vector<std::string>> rows(num_events, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_events; i++)
	{
		rows[i][0] = stamp2string(l2[i].Tim_s, l2[i].Tim_ms);
		rows[i][1] = ss_double(l2[i].Tim_s+l2[i].Tim_ms/1000. + 1609459200);
		rows[i][2] = ss_double(l2[i].Star_Z_Solar_angle);
		rows[i][3] = ss_double(l2[i].Moon_J2000_q1);
		rows[i][4] = ss_double(l2[i].Moon_J2000_q2);
		rows[i][5] = ss_double(l2[i].Moon_J2000_q3);
		rows[i][6] = ss_int(l2[i].Moon_J2000_Star_X);
		rows[i][7] = ss_int(l2[i].Moon_J2000_Star_Y);
		rows[i][8] = ss_int(l2[i].Moon_J2000_Star_Z);
		rows[i][9] = ss_double(l2[i].fly_time);
		rows[i][10] = ss_double(l2[i].x1);
		rows[i][11] = ss_double(l2[i].y1);
		rows[i][12] = ss_double(l2[i].x2);
		rows[i][13] = ss_double(l2[i].y2);
		rows[i][14] = ss_double(l2[i].energy/1000.);
		rows[i][15] = (l2[i].spec == 0) ? ("H") : ("O");
		rows[i][16] = ss_double(l2[i].ep_mv);
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}
int OutputL0_csv_p(data_ary<GenaL0> ary, const std::string &fileName)
{
	GenaL0 *l0 = ary.data;
	size_t num_packets = ary.count;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);
	int columns = 17;
	// 初始化二维vector，外层vector大小为num * 62，代表总共有num * 62行数据
	std::vector<std::vector<std::string>> rows(num_packets * 62, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_packets * 62; i++)
	{
		// 计算当前属于第几组数据（每组对应62行）
		size_t packet_idx = i / 62;
		// 计算在当前组内是第几行
		int row_in_packet = i % 62;
		// BJ time
        rows[i][0] = stamp2string(l0[packet_idx].Tim_s, l0[packet_idx].Tim_ms);
        rows[i][1] = ss_double(l0[packet_idx].Tim_s + l0[packet_idx].Tim_ms/1000.);
		rows[i][2] = ss_double(l0[packet_idx].evts[row_in_packet].internal_time);

		for (int j = 0; j < 4; j++)
		{
			rows[i][3 + j] = ss_int(l0[packet_idx].evts[row_in_packet].tof[j]);
		}
		for (int j = 0; j < 2; j++)
		{
			rows[i][7 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c1[j]);
			rows[i][9 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c3[j]);
			rows[i][11 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c5[j]);
			rows[i][13 + j] = ss_int(l0[packet_idx].evts[row_in_packet].c7[j]);
		}
		rows[i][15] = ss_int(l0[packet_idx].evts[row_in_packet].ep);
		rows[i][16] = ss_int(l0[packet_idx].evts[row_in_packet].isfill);
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}


int OutputL1_csv_p(data_ary<GenaL1> ary, const std::string &fileName)
{
	GenaL1 *l1 = ary.data;
	size_t num_events = ary.count;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);

	int columns = 15;
	std::vector<std::vector<std::string>> rows(num_events, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_events; i++)
	{
		rows[i][0] = stamp2string(l1[i].Tim_s, l1[i].Tim_ms);
		rows[i][1] = ss_double(l1[i].Tim_s+l1[i].Tim_ms/1000.+ 1609459200);
		for (int j = 0; j < 4; ++j)
		{
			rows[i][2 + j] = ss_double(l1[i].fly_time[j]);
		}
	
		rows[i][6] = ss_double(l1[i].tx1[0]);
		rows[i][7] = ss_double(l1[i].tx2[0]);
		
        rows[i][8] = ss_double(l1[i].ty1[0]);
		rows[i][9] = ss_double(l1[i].ty2[0]);

		rows[i][10] = ss_double(l1[i].tx1[1]);
		rows[i][11] = ss_double(l1[i].tx2[1]);
		
        rows[i][12] = ss_double(l1[i].ty1[1]);
		rows[i][13] = ss_double(l1[i].ty2[1]);
		
        rows[i][14] = ss_double(l1[i].ep_mv);
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}

int OutputL2_csv_p(data_ary<GenaL2> ary, const std::string &fileName)
{
	GenaL2 *l2 = ary.data;
	size_t num_events = ary.count;
	std::ofstream stream(fileName, std::ios::app);
	Writer<delimiter<','>> writer(stream);

	int columns = 10;
	std::vector<std::vector<std::string>> rows(num_events, std::vector<std::string>(columns));
	for (size_t i = 0; i < num_events; i++)
	{

		rows[i][0] = stamp2string(l2[i].Tim_s, l2[i].Tim_ms);
		rows[i][1] = ss_double(l2[i].Tim_s+l2[i].Tim_ms/1000.+ 1609459200);
		rows[i][2] = ss_double(l2[i].fly_time);
		rows[i][3] = ss_double(l2[i].x1);
		rows[i][4] = ss_double(l2[i].y1);
		rows[i][5] = ss_double(l2[i].x2);
		rows[i][6] = ss_double(l2[i].y2);
		rows[i][7] = ss_double(l2[i].energy/1000.);
		rows[i][8] = (l2[i].spec == 0) ? ("H") : ("O");
		rows[i][9] = ss_double(l2[i].ep_mv);
	}
	writer.write_rows(rows);
	stream.close();
	return 0;
}

// 函数用于判断字符串是否全由空白字符组成
bool isAllWhitespace(const char* str)
{
    while (*str!= '\0')
    {
        if (!isspace(*str))  // 如果当前字符不是空白字符（空格、制表符、换行等）
        {
            return false;
        }
        str++;
    }
    return true;
}
