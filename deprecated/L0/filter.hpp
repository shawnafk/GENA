#ifndef FILTER_HPP
#define FILTER_HPP
#include <utility>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include "structures.hpp"
using namespace std;

// 模板函数，用于根据给定的筛选规则筛选数据，元素类型和数组类型都使用模板参数
template<typename T, typename Func>
inline data_ary<T> filter_rules(const data_ary<T>& arrayInfo, Func rule) {
    std::vector<int> idx;
    // 根据传入的筛选规则函数rule来获取符合条件的索引，遍历数组
    for (size_t i = 0; i < arrayInfo.count; ++i) {
        if (rule(arrayInfo.data[i])) {
            idx.push_back(i);
        }
    }

    // 开辟一段内存空间用于存储筛选后的结果数据，大小为符合条件元素个数
    T* result = static_cast<T*>(std::malloc(idx.size() * sizeof(T)));
    if (!result) {
        // 如果内存分配失败，进行错误处理，这里简单返回默认初始化的data_ary，实际可根据需求完善
        return {nullptr, 0};
    }

    size_t index = 0;
    // 根据索引从原始数据中提取出对应的元素，拷贝到新开辟的内存空间中
    for (int i : idx) {
        new (&result[index]) T(arrayInfo.data[i]);  // 使用placement new来在指定内存位置构造对象
        ++index;
    }

    // 返回包含筛选后数据指针和元素数量的结构体
    return {result, idx.size()};
}
bool rule_prob_a(const GenaL1 &l1);
bool rule_prob_b(const GenaL1 &l1);
bool rule_L1A(const GenaL1 &l1);
bool rule_L1B(const GenaL1 &l1);
bool rule_L1C(const GenaL1 &l1);
bool rule_L1D(const GenaL1 &l1);
bool rule_L2A(const GenaL2 &l2);
bool rule_L2B(const GenaL2 &l2);
data_ary<GenaL0> raw_2_L0( Packet *packet, unsigned long long num_patch, unsigned long long max_packets);
data_ary<GenaL1> L0_2_L1(data_ary<GenaL0>);
data_ary<GenaL2> L1_2_L2(data_ary<GenaL1>);

//pair<data_ary<GenaParameters_A>,data_ary<GenaParameters_B>> para_from_raw(Packet *packet, unsigned long long num_patch);
data_ary<GenaParameters> para_from_raw(Packet *packet, unsigned long long num_patch);

#endif
