把comparator的构造放在了tools.cpp中
但是把context的构造放在了generate_comparator外面
因为会涉及动态内存分配的问题

秘密共享拆分的时候
mod的是数据，而不是数据范围
保证拆分值都是整数

计算方差遗留问题
可以考虑ciphertext packing
也可以就单独加密{01011100...}
先用最基本的逐比特加密把！

其实理论上来说
cloud2上存放的kd_tree中的划分信息N都是要经过shuffle的
但是本人，为了图方便，现在暂且先略过
后期可能就要加上
cloud2的N加密后发给cloud1要经过shuffle还原才能继续操作
蒜了蒜了，到时候做的好再加把
万一糊了，还能少留点泪水