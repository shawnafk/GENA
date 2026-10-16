#ifndef IO_HPP
#define IO_HPP
std::string ss_int(int num);
std::string ss_double(double num);
int Init_csv(const std::string& fileName, int level);
int Init_csv_p(const std::string& fileName, int level);
int Init_csv_parameter(const std::string &fileName);
int Output_parameter_csv(data_ary<GenaParameters> ary, const std::string &fileName);
int OutputL0_csv(data_ary<GenaL0> l0, const std::string& fileName);
int OutputL1_csv(data_ary<GenaL1> l1, const std::string& fileName);
int OutputL2_csv(data_ary<GenaL2> l2, const std::string& fileName);
int OutputL0_csv_p(data_ary<GenaL0> l0, const std::string& fileName);
int OutputL1_csv_p(data_ary<GenaL1> l1, const std::string& fileName);
int OutputL2_csv_p(data_ary<GenaL2> l2, const std::string& fileName);
bool isAllWhitespace(const char* str);
#endif
