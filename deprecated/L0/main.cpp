#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <algorithm>
// #include <omp.h>
#include "structures.hpp"
#include "io.hpp"
#include "packets.hpp"
#include "filter.hpp"
#include <Python.h>

#define PATCH_SIZE 30000
using namespace std;
#include <filesystem>
#include <dirent.h>
namespace fs = std::filesystem;

double A_X_LEFT_NS;
double A_X_RIGHT_NS;
double A_Y_LEFT_NS;
double A_Y_RIGHT_NS;

double B_X_LEFT_NS;
double B_X_RIGHT_NS;
double B_Y_LEFT_NS;
double B_Y_RIGHT_NS;

double A_NS_XLIM;
double A_NS_YLIM;
double A_X_MAX_NS;
double A_Y_MAX_NS;

double A_X_CENTER_NS;
double A_Y_CENTER_NS;

double A_X_COF;
double A_Y_COF;

double B_X_MAX_NS;
double B_Y_MAX_NS;
double B_NS_XLIM;
double B_NS_YLIM;

double B_X_CENTER_NS;
double B_Y_CENTER_NS;

double B_X_COF;
double B_Y_COF;
// 1: orbit, 2: ground
int SN_ID;

int main(int argc, char *argv[])
{
	clock_t begin = clock();

	Py_Initialize();
	if (!Py_IsInitialized())
	{
		std::cerr << "Failed to initialize Python interpreter" << std::endl;
		return 1;
	}

	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('./')");

	// argv[1] is config file name
	PyObject *pModule = PyImport_ImportModule(argv[1]);
	if (pModule == nullptr)
	{
		std::cerr << "Failed to load config" << std::endl;
		Py_Finalize();
		return 1;
	}
	// GLOBAL VARIABLES
	SN_ID = PyLong_AsLong(PyObject_GetAttrString(pModule, "SN_ID"));

	A_X_LEFT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "A_X_LEFT_NS"));
	A_Y_LEFT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "A_Y_LEFT_NS"));
	A_X_RIGHT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "A_X_RIGHT_NS"));
	A_Y_RIGHT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "A_Y_RIGHT_NS"));

	B_X_LEFT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "B_X_LEFT_NS"));
	B_Y_LEFT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "B_Y_LEFT_NS"));
	B_X_RIGHT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "B_X_RIGHT_NS"));
	B_Y_RIGHT_NS = PyFloat_AsDouble(PyObject_GetAttrString(pModule, "B_Y_RIGHT_NS"));

	A_X_MAX_NS = (abs(A_X_RIGHT_NS) + abs(A_X_LEFT_NS));
	A_Y_MAX_NS = (abs(A_Y_RIGHT_NS) + abs(A_Y_LEFT_NS));
	// need divided by 2 and yes
	A_NS_XLIM = A_X_MAX_NS / 2;
	A_NS_YLIM = A_Y_MAX_NS / 2;
	A_X_CENTER_NS = (A_X_LEFT_NS + A_X_RIGHT_NS);
	A_Y_CENTER_NS = (A_Y_LEFT_NS + A_Y_RIGHT_NS);
	A_X_COF = X_MCP / A_X_MAX_NS;
	A_Y_COF = Y_MCP / A_Y_MAX_NS;

	B_X_MAX_NS = abs(B_X_RIGHT_NS) + abs(B_X_LEFT_NS);
	B_Y_MAX_NS = abs(B_Y_RIGHT_NS) + abs(B_Y_LEFT_NS);
	B_NS_XLIM = B_X_MAX_NS / 2;
	B_NS_YLIM = B_Y_MAX_NS / 2;
	// should divided by 2 but not yet
	B_X_CENTER_NS = (B_X_LEFT_NS + B_X_RIGHT_NS);
	B_Y_CENTER_NS = (B_Y_LEFT_NS + B_Y_RIGHT_NS);
	B_X_COF = X_MCP / B_X_MAX_NS;
	B_Y_COF = Y_MCP / B_Y_MAX_NS;

	//---------------
	// load parameters from python file

	std::vector<std::string> c_InputList;
	PyObject *InputList = PyObject_GetAttrString(pModule, "InputList");
	Py_ssize_t listSize = PyList_Size(InputList);
	for (Py_ssize_t i = 0; i < listSize; ++i)
	{
		PyObject *item = PyList_GetItem(InputList, i);
		PyObject *tempBytes = PyUnicode_AsEncodedString(item, "utf-8", "strict");
		std::string value = PyBytes_AsString(tempBytes);
		c_InputList.push_back(value);
		Py_DECREF(tempBytes);
		Py_XDECREF(item);
	}

	Py_XDECREF(InputList);

	int partial = PyLong_AsLong(PyObject_GetAttrString(pModule, "PartialOutput"));
	int outputLevel0 = PyLong_AsLong(PyObject_GetAttrString(pModule, "L0"));
	int outputLevel1 = PyLong_AsLong(PyObject_GetAttrString(pModule, "L1"));
	int outputLevel1a = PyLong_AsLong(PyObject_GetAttrString(pModule, "L1A"));
	int outputLevel1b = PyLong_AsLong(PyObject_GetAttrString(pModule, "L1B"));
	int outputLevel1c = PyLong_AsLong(PyObject_GetAttrString(pModule, "L1C"));
	int outputLevel1d = PyLong_AsLong(PyObject_GetAttrString(pModule, "L1D"));
	int outputLevel2 = PyLong_AsLong(PyObject_GetAttrString(pModule, "L2"));
	int outputLevel2a = PyLong_AsLong(PyObject_GetAttrString(pModule, "L2A"));
	int outputLevel2b = PyLong_AsLong(PyObject_GetAttrString(pModule, "L2B"));

	// PyObject* tempBytes = PyUnicode_AsEncodedString(PyObject_GetAttrString(pModule, "OutputFolder"), "utf-8", "strict");
	// std::string outpath = PyBytes_AsString(tempBytes);

	int parse_data = PyLong_AsLong(PyObject_GetAttrString(pModule, "Parse_data"));
	int parse_parameter = PyLong_AsLong(PyObject_GetAttrString(pModule, "Parse_parameter"));
	// 关闭 Python 解释器

	Py_DECREF(pModule);
	Py_Finalize();

	//---------------------------------------------
	//	if (!fs::exists(outpath))
	//	{
	//		try
	//		{
	//			fs::create_directories(outpath);
	//			std::cout << "create path success" << std::endl;
	//		}
	//		catch (const fs::filesystem_error &e)
	//		{
	//			std::cerr << "create path failed: " << e.what() << std::endl;
	//			return 1;
	//		}
	//	}
	//	else
	//	{
	//		std::cout << "path already exists" << std::endl;
	//	}

	clock_t time_start, time_end;
	double elapsed_time;
	// loop raw files
	for (int j = 0; j < c_InputList.size(); j++)
	{
		// every file in the input list has its own output name
		string full_outpath = c_InputList[j].data() + std::string("_");
		cout << full_outpath << endl;

		time_start = clock();
		cout << "Processing file " << c_InputList[j] << endl;
		FILE *rawfile = fopen(c_InputList[j].data(), "rb");
		if (rawfile == NULL)
		{
			perror("can not open file");
			// continue;
		}
		else
		{
			long long fstart, fend;
			if (0)
			{
				fstart = 2064L * 1000L;
				fend = fstart + 2064L * 2000L;
			}
			else
			{
				fstart = 0;
				fseek(rawfile, 0, SEEK_END);
				fend = ftell(rawfile);
			}
			if (fend == -1)
			{
				std::cerr << "open file failed!" << std::endl;
				return -1;
			}
			std::cout << "file length: " << fend - fstart << std::endl;
			long buffer_size = fend - fstart;
			uint8_t *raw_buffer = (uint8_t *)malloc(buffer_size * sizeof(uint8_t));
			if (raw_buffer == NULL)
			{
				std::cerr << "alloc raw_buffer failed!" << std::endl;
				return 1;
			}
			fseek(rawfile, 0, SEEK_SET);
			if (fread(raw_buffer, sizeof(uint8_t), buffer_size, rawfile) == 0)
			{
				std::cerr << "read to buffer failed" << std::endl;
			}
			fclose(rawfile);

			time_end = clock();
			elapsed_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
			std::cout << "read raw to buffer cost: " << elapsed_time << " s" << std::endl;

			time_start = clock();
			std::vector<std::int64_t> all_header, selected_header;
			long long len_range[2] = {fstart, fend};
			all_header = get_packet_header_addr(raw_buffer, len_range);
			time_end = clock();
			elapsed_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
			std::cout << "find header addr cost: " << elapsed_time << " s" << std::endl;
            
			time_start = clock();
            selected_header = find_2064(all_header, len_range[1] - len_range[0]);
            long long num_packets = selected_header.size();
			std::cout << "found packets num: " << num_packets << std::endl;
			Packet *packet = (Packet *)malloc(num_packets * sizeof(Packet));
			// std::cout << "sizeof Packet: " << sizeof(Packet) << std::endl;
			{
				// from file to complete packets
				for (long i = 0; i < num_packets; i++)
				{
					std::memcpy(&packet[i], raw_buffer + selected_header[i], sizeof(Packet));
				}
			}
			free(raw_buffer);
			time_end = clock();
			elapsed_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
			std::cout << "select and load all packets cost: " << elapsed_time << " s" << std::endl;

			// read engineering parameters from all packets

			if (parse_parameter == 1)
			{
				Init_csv_parameter(full_outpath + "parameter.csv");
				
				time_start = clock();
				data_ary<GenaParameters> ary_para;
				ary_para = para_from_raw(packet, num_packets);
				time_end = clock();
				elapsed_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
				std::cout << "parse enigeering parameters cost: " << elapsed_time << " s" << std::endl;
				
				Output_parameter_csv(ary_para, full_outpath + "parameter.csv");
			}

			if (parse_data == 1)
			{

				// if parse data then init csv

				if (partial)
				{
					std::cout << "partial output" << std::endl;
					if (outputLevel0)
					{
						// L0, measured raw (pass)
						Init_csv_p(full_outpath + "L0.csv", 0);
						std::cout << "output 0 level " << std::endl;
					}
					if (outputLevel1)
					{
						// L1, to ns, and mV (pass)
						Init_csv_p(full_outpath + "L1_a.csv", 1);
						Init_csv_p(full_outpath + "L1_b.csv", 1);
						std::cout << "output 1 level " << std::endl;
					}
					if (outputLevel1a)
					{
						// 1A select anode not all zero (pass)
						Init_csv_p(full_outpath + "L1A_a.csv", 1);
						Init_csv_p(full_outpath + "L1A_b.csv", 1);
						std::cout << "output 1A level " << std::endl;
					}
					if (outputLevel1b)
					{
						// 1B (used to obtain L2 data) 4 start signal with ns limited; 4 stop all exists. (pass)
						Init_csv_p(full_outpath + "L1B_a.csv", 1);
						Init_csv_p(full_outpath + "L1B_b.csv", 1);
						std::cout << "output 1B level " << std::endl;
					}
					if (outputLevel1c)
					{
						// 1B (used to obtain L2 data) 4 start signal with ns limited; 4 stop all exists. (pass)
						Init_csv_p(full_outpath + "L1C_a.csv", 1);
						Init_csv_p(full_outpath + "L1C_b.csv", 1);
						std::cout << "output 1C level " << std::endl;
					}
					if (outputLevel1d)
					{
						// 1B (used to obtain L2 data) 4 start signal with ns limited; 4 stop all exists. (pass)
						Init_csv_p(full_outpath + "L1D_a.csv", 1);
						Init_csv_p(full_outpath + "L1D_b.csv", 1);
						std::cout << "output 1D level " << std::endl;
					}
					if (outputLevel2)
					{
						// L2, direct convert ns to mm, with no other restriction
						Init_csv_p(full_outpath + "L2_a.csv", 2);
						Init_csv_p(full_outpath + "L2_b.csv", 2);
						std::cout << "output 2 level " << std::endl;
					}
					if (outputLevel2a)
					{
						// L2, direct convert ns to mm, with no other restriction
						Init_csv_p(full_outpath + "L2A_a.csv", 2);
						Init_csv_p(full_outpath + "L2A_b.csv", 2);
						std::cout << "output 2A level " << std::endl;
					}
					if (outputLevel2b)
					{
						// L2, direct convert ns to mm, with no other restriction
						Init_csv_p(full_outpath + "L2B_a.csv", 2);
						Init_csv_p(full_outpath + "L2B_b.csv", 2);
						std::cout << "output 2B level " << std::endl;
					}
				}
				else
				{
					std::cout << "full output" << std::endl;
					if (outputLevel0)
					{
						// L0, measured raw (pass)
						Init_csv(full_outpath + "L0.csv", 0);
						std::cout << "output 0 level " << std::endl;
					}
					if (outputLevel1)
					{
						// L1, to ns, and mV (pass)
						Init_csv(full_outpath + "L1_a.csv", 1);
						Init_csv(full_outpath + "L1_b.csv", 1);
						std::cout << "output 1 level " << std::endl;
					}
					if (outputLevel1a)
					{
						// 1A select anode not all zero (pass)
						Init_csv(full_outpath + "L1A_a.csv", 1);
						Init_csv(full_outpath + "L1A_b.csv", 1);
						std::cout << "output 1A level " << std::endl;
					}
					if (outputLevel1b)
					{
						// 1B (used to obtain L2 data) 4 start signal with ns limited; 4 stop all exists. (pass)
						Init_csv(full_outpath + "L1B_a.csv", 1);
						Init_csv(full_outpath + "L1B_b.csv", 1);
						std::cout << "output 1B level " << std::endl;
					}
					if (outputLevel1c)
					{
						// 1B (used to obtain L2 data) 4 start signal with ns limited; 4 stop all exists. (pass)
						Init_csv(full_outpath + "L1C_a.csv", 1);
						Init_csv(full_outpath + "L1C_b.csv", 1);
						std::cout << "output 1C level " << std::endl;
					}
					if (outputLevel1d)
					{
						// 1B (used to obtain L2 data) 4 start signal with ns limited; 4 stop all exists. (pass)
						Init_csv(full_outpath + "L1D_a.csv", 1);
						Init_csv(full_outpath + "L1D_b.csv", 1);
						std::cout << "output 1D level " << std::endl;
					}
					if (outputLevel2)
					{
						// L2, direct convert ns to mm, with no other restriction
						Init_csv(full_outpath + "L2_a.csv", 2);
						Init_csv(full_outpath + "L2_b.csv", 2);
						std::cout << "output 2 level " << std::endl;
					}
					if (outputLevel2a)
					{
						// L2, direct convert ns to mm, with no other restriction
						Init_csv(full_outpath + "L2A_a.csv", 2);
						Init_csv(full_outpath + "L2A_b.csv", 2);
						std::cout << "output 2A level " << std::endl;
					}
					if (outputLevel2b)
					{
						// L2, direct convert ns to mm, with no other restriction
						Init_csv(full_outpath + "L2B_a.csv", 2);
						Init_csv(full_outpath + "L2B_b.csv", 2);
						std::cout << "output 2B level " << std::endl;
					}
				}

				// parse data in patches, save memory!
				time_start = clock();
				unsigned long remaining_packets = num_packets;
				unsigned long num_patch = 0;
				unsigned long cur_packet = 0;
				unsigned long patch_id = 0;
				while (remaining_packets > 0)
				{
					// std::cout << "processing " << (1 - double(remaining_packets) / double(num_packets)) * 100 << "%" << endl;
					num_patch = remaining_packets < PATCH_SIZE ? remaining_packets : PATCH_SIZE;

					data_ary<GenaL0> ary_L0 = raw_2_L0(&packet[cur_packet], num_patch, num_packets);
					cur_packet = cur_packet + num_patch;
					if (outputLevel0)
					{
						if (partial)
						{
							// L0, measured raw (pass)
							OutputL0_csv_p(ary_L0, full_outpath + "L0.csv");
						}
						else
						{
							OutputL0_csv(ary_L0, full_outpath + "L0.csv");
						}
					}

					// ------------------------------------------------------------- FILTER and OUTPUT
					time_start = clock();
					data_ary<GenaL1> ary_L1 = L0_2_L1(ary_L0);
					free(ary_L0.data);

					//----------------------------------------------selection ---------------------------------------
					// all start from L1 A and B, independtly
					data_ary<GenaL1> ary_L1_a = filter_rules(ary_L1, rule_prob_a);
					// std::cout << "A num is " << ary_L1_a.count << std::endl;
					data_ary<GenaL1> ary_L1_b = filter_rules(ary_L1, rule_prob_b);
					// std::cout << "B num is " << ary_L1_b.count << std::endl;
					free(ary_L1.data);

					if (outputLevel1)
					{
						// L1, to ns, and mV (pass)
						if (partial)
						{
							OutputL1_csv_p(ary_L1_a, full_outpath + "L1_a.csv");
							OutputL1_csv_p(ary_L1_b, full_outpath + "L1_b.csv");
						}
						else
						{
							OutputL1_csv(ary_L1_a, full_outpath + "L1_a.csv");
							OutputL1_csv(ary_L1_b, full_outpath + "L1_b.csv");
						}
					}

					if (outputLevel1a)
					{
						data_ary<GenaL1> ary_L1A_a = filter_rules(ary_L1_a, rule_L1A);
						data_ary<GenaL1> ary_L1A_b = filter_rules(ary_L1_b, rule_L1A);
						if (partial)
						{
							OutputL1_csv_p(ary_L1A_a, full_outpath + "L1A_a.csv");
							OutputL1_csv_p(ary_L1A_b, full_outpath + "L1A_b.csv");
						}
						else
						{
							OutputL1_csv(ary_L1A_a, full_outpath + "L1A_a.csv");
							OutputL1_csv(ary_L1A_b, full_outpath + "L1A_b.csv");
						}
						free(ary_L1A_a.data);
						free(ary_L1A_b.data);
					}

					if (outputLevel1b)
					{
						data_ary<GenaL1> ary_L1B_a = filter_rules(ary_L1_a, rule_L1B);
						data_ary<GenaL1> ary_L1B_b = filter_rules(ary_L1_b, rule_L1B);
						// std::cout << "count 1B_a = " << ary_L1B_a.count << std::endl;
						// std::cout << "count 1B_b = " << ary_L1B_b.count << std::endl;
						if (partial)
						{
							OutputL1_csv_p(ary_L1B_a, full_outpath + "L1B_a.csv");
							OutputL1_csv_p(ary_L1B_b, full_outpath + "L1B_b.csv");
						}
						else
						{
							OutputL1_csv(ary_L1B_a, full_outpath + "L1B_a.csv");
							OutputL1_csv(ary_L1B_b, full_outpath + "L1B_b.csv");
						}
						free(ary_L1B_a.data);
						free(ary_L1B_b.data);
					}

					if (outputLevel1c)
					{
						data_ary<GenaL1> ary_L1C_a = filter_rules(ary_L1_a, rule_L1C);
						data_ary<GenaL1> ary_L1C_b = filter_rules(ary_L1_b, rule_L1C);
						if (partial)
						{
							OutputL1_csv_p(ary_L1C_a, full_outpath + "L1C_a.csv");
							OutputL1_csv_p(ary_L1C_b, full_outpath + "L1C_b.csv");
						}
						else
						{
							OutputL1_csv(ary_L1C_a, full_outpath + "L1C_a.csv");
							OutputL1_csv(ary_L1C_b, full_outpath + "L1C_b.csv");
						}
						free(ary_L1C_a.data);
						free(ary_L1C_b.data);
					}

					data_ary<GenaL1> ary_L1D_a = filter_rules(ary_L1_a, rule_L1D);
					data_ary<GenaL1> ary_L1D_b = filter_rules(ary_L1_b, rule_L1D);

					if (outputLevel1d)
					{
						if (partial)
						{
							OutputL1_csv_p(ary_L1D_a, full_outpath + "L1D_a.csv");
							OutputL1_csv_p(ary_L1D_b, full_outpath + "L1D_b.csv");
						}
						else
						{
							OutputL1_csv(ary_L1D_a, full_outpath + "L1D_a.csv");
							OutputL1_csv(ary_L1D_b, full_outpath + "L1D_b.csv");
						}
					}

					data_ary<GenaL2> ary_L2_a = L1_2_L2(ary_L1D_a);
					data_ary<GenaL2> ary_L2_b = L1_2_L2(ary_L1D_b);
					free(ary_L1D_a.data);
					free(ary_L1D_b.data);
					free(ary_L1_a.data);
					free(ary_L1_b.data);

					// to the next level
					if (outputLevel2)
					{

						if (partial)
						{
							OutputL2_csv_p(ary_L2_a, full_outpath + "L2_a.csv");
							OutputL2_csv_p(ary_L2_b, full_outpath + "L2_b.csv");
						}
						else
						{
							OutputL2_csv(ary_L2_a, full_outpath + "L2_a.csv");
							OutputL2_csv(ary_L2_b, full_outpath + "L2_b.csv");
						}
					}

					if (outputLevel2a)
					{
						// select from L1 and convert to L2
						data_ary<GenaL2> ary_L2A_a = filter_rules(ary_L2_a, rule_L2A);
						data_ary<GenaL2> ary_L2A_b = filter_rules(ary_L2_b, rule_L2A);
						if (partial)
						{
							OutputL2_csv_p(ary_L2A_a, full_outpath + "L2A_a.csv");
							OutputL2_csv_p(ary_L2A_b, full_outpath + "L2A_b.csv");
						}
						else
						{
							OutputL2_csv(ary_L2A_a, full_outpath + "L2A_a.csv");
							OutputL2_csv(ary_L2A_b, full_outpath + "L2A_b.csv");
						}
						free(ary_L2A_a.data);
						free(ary_L2A_b.data);
					}

					if (outputLevel2b)
					{
						// select from L1 and convert to L2
						data_ary<GenaL2> ary_L2B_a = filter_rules(ary_L2_a, rule_L2B);
						data_ary<GenaL2> ary_L2B_b = filter_rules(ary_L2_b, rule_L2B);
						if (partial)
						{
							OutputL2_csv_p(ary_L2B_a, full_outpath + "L2B_a.csv");
							OutputL2_csv_p(ary_L2B_b, full_outpath + "L2B_b.csv");
						}
						else
						{
							OutputL2_csv(ary_L2B_a, full_outpath + "L2B_a.csv");
							OutputL2_csv(ary_L2B_b, full_outpath + "L2B_b.csv");
						}
						free(ary_L2B_a.data);
						free(ary_L2B_b.data);
					}

					free(ary_L2_a.data);
					free(ary_L2_b.data);
					remaining_packets = remaining_packets - num_patch;
					patch_id = patch_id + 1;
					//time_end = clock();
					//elapsed_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
					//cout << "Filter and IO" << " cost:" << elapsed_time << " s" << endl;
				}
				 cout << "processing 100%" << endl;
				free(packet);
				time_end = clock();
				elapsed_time = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
				std::cout << "parse packets cost: " << elapsed_time << " s" << std::endl;
			}
		}
	}
	clock_t over = clock();
	elapsed_time = static_cast<double>(over - begin) / CLOCKS_PER_SEC;
	cout << "total time cost:" << elapsed_time << " s" << endl;
	return 0;
}
