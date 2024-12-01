

# A fast and elitist multiobjective genetic algorithm: NSGA-II

Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. A. M. T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE transactions on evolutionary computation, 6(2), 182-197. 

> Homework of Optimization-Method class

## 算法流程

1. 生成初始种群作为亲本（N）
2. 通过遗传算法对亲本种群进行选择、交叉、变异，得到子代种群（2N）
3. 将亲代和子代种群合并，对合并后得种群进行非支配排序和拥挤度计算
4. 通过精英选择策略，从合并的种群中选择N个个体作为新的亲本
5. 循环迭代直到设定的终止条件

## 对应函数

- 主函数（main）
- 初始化种群（initialize_variables）
- 快速非支配排序和拥挤度计算（non_domination_sort_mod）
- 竞标赛选择代码（tournament_selection）
- 交叉 变异代码（genetic_operator）
- 生成新的种群（精英策略）（replace_chromosome）
- 测试函数（evaluate_objective）

> 
