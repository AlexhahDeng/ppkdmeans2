#include "tools.h"

//! 自定义工具函数
unsigned long p = 7;            // plaintext modulus 明文模数
unsigned long d = 2;            // field extension degree 域扩展次数
unsigned long m = 300;          // cyclotomic polynomial 分圆多项式次数
unsigned long nb_primes = 170;  // modulus chain中密文素数比特的数量
unsigned long c = 3;            // 密钥切换矩阵的列数


auto context = ContextBuilder<BGV>()// 初始化上下文
    .m(m)// 分圆环的规模
    .p(p)// 明文素数模数
    .r(1)// 默认的什么东西
    .bits(nb_primes)// 模数链中密文素数比特的数量
    .c(c)// 密钥转换矩阵的列数
    .scale(6)// 不知道啥东西
    .build();

//! 密文比较部分的tools
void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots)
{
	decomp.clear();
	decomp.resize(nslots, 0);
	int power = intlog(base, input) + 1;
	if (power > nslots)
	{
		cout << "Input character is too big to be converted" << endl;
		exit(1);
	}
	unsigned long rest = input;
	unsigned long coeff;

	int i = 0;
	while (i < power)
	{
		coeff = rest % base;
		decomp[i] = coeff;
		rest = (rest - coeff) / base;
		i++;
	}
}

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
void simplePolyEval(Ctxt& ret, const NTL::ZZX& poly, DynamicCtxtPowers& babyStep)
{
	ret.clear();
	if (deg(poly) < 0) return;       // the zero polynomial always returns zero

	//OLD: assert(deg(poly)<=babyStep.size()); // ensure that we have enough powers
	helib::assertTrue(deg(poly) <= babyStep.size(), "BabyStep has not enough powers (required more than deg(poly))");

	NTL::ZZ coef;
	NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
	for (long i = 1; i <= deg(poly); i++) {
		rem(coef, coeff(poly, i), p);
		if (coef > p / 2) coef -= p;

		Ctxt tmp = babyStep.getPower(i); // X^i
		tmp.multByConstant(coef);        // f_i X^i
		ret += tmp;
	}
	// Add the free term
	rem(coef, ConstTerm(poly), p);
	if (coef > p / 2) coef -= p;
	ret.addConstant(coef);
	//  if (verbose) checkPolyEval(ret, babyStep[0], poly);
}

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
void PatersonStockmeyer(Ctxt& ret, const NTL::ZZX& poly, long k, long t, long delta,
	DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep)
{
	if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
		simplePolyEval(ret, poly, babyStep);
		return;
	}
	NTL::ZZX r = trunc(poly, k * t);      // degree <= k*2^e-1
	NTL::ZZX q = RightShift(poly, k * t); // degree == k(2^e-1) +delta

	const NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
	const NTL::ZZ& coef = coeff(r, deg(q));
	SetCoeff(r, deg(q), coef - 1);  // r' = r - X^{deg(q)}

	NTL::ZZX c, s;
	DivRem(c, s, r, q); // r' = c*q + s
	// deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

	helib::assertTrue(deg(s) < deg(q), "Degree of s is not less than degree of q");
	helib::assertTrue(IsZero(c) || deg(c) < k - delta, "Nonzero c has not degree smaller than k - delta");
	SetCoeff(s, deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

	// reduce the coefficients modulo p
	for (long i = 0; i <= deg(c); i++) rem(c[i], c[i], p);
	c.normalize();
	for (long i = 0; i <= deg(s); i++) rem(s[i], s[i], p);
	s.normalize();

	// Evaluate recursively poly = (c+X^{kt})*q + s'
	PatersonStockmeyer(ret, q, k, t / 2, delta, babyStep, giantStep);

	Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
	simplePolyEval(tmp, c, babyStep);
	tmp += giantStep.getPower(t);
	ret.multiplyBy(tmp);

	PatersonStockmeyer(tmp, s, k, t / 2, delta, babyStep, giantStep);
	ret += tmp;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
void degPowerOfTwo(Ctxt& ret, const NTL::ZZX& poly, long k,
	DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep)
{
	if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
		simplePolyEval(ret, poly, babyStep);
		return;
	}
	long n = divc(deg(poly), k);        // We assume n=2^e or n=2^e -1
	n = 1L << NTL::NextPowerOfTwo(n); // round up to n=2^e
	NTL::ZZX r = trunc(poly, (n - 1) * k);      // degree <= k(2^e-1)-1
	NTL::ZZX q = RightShift(poly, (n - 1) * k); // 0 < degree < 2k
	SetCoeff(r, (n - 1) * k);              // monic, degree == k(2^e-1)
	q -= 1;

	PatersonStockmeyer(ret, r, k, n / 2, 0, babyStep, giantStep);

	Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
	simplePolyEval(tmp, q, babyStep); // evaluate q

	// multiply by X^{k(n-1)} with minimum depth
	for (long i = 1; i < n; i *= 2) {
		tmp.multiplyBy(giantStep.getPower(i));
	}
	ret += tmp;
}

void recursivePolyEval(Ctxt& ret, const NTL::ZZX& poly, long k,
	DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep)
{
	if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
		simplePolyEval(ret, poly, babyStep);
		return;
	}

	long delta = deg(poly) % k; // deg(poly) mod k
	long n = divc(deg(poly), k); // ceil( deg(poly)/k )
	long t = 1L << (NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

	// Special case for deg(poly) = k * 2^e +delta
	if (n == t) {
		degPowerOfTwo(ret, poly, k, babyStep, giantStep);
		return;
	}

	// When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
	if (n == t - 1 && delta == 0) {
		PatersonStockmeyer(ret, poly, k, t / 2, delta, babyStep, giantStep);
		return;
	}

	t = t / 2;

	// In any other case we have kt < deg(poly) < k(2t-1). We then set 
	// u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
	// and recurse on poly = (q-1)*X^u + (X^u+r)

	long u = deg(poly) - k * (t - 1);
	NTL::ZZX r = trunc(poly, u);      // degree <= u-1
	NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
	q -= 1;
	SetCoeff(r, u);              // degree == u

	PatersonStockmeyer(ret, q, k, t / 2, 0, babyStep, giantStep);

	Ctxt tmp = giantStep.getPower(u / k);
	if (delta != 0) { // if u is not divisible by k then compute it
		tmp.multiplyBy(babyStep.getPower(delta));
	}
	ret.multiplyBy(tmp);

	recursivePolyEval(tmp, r, k, babyStep, giantStep);
	ret += tmp;
}

//! tools for data processing
void read_data() {
	// 29 attributes, 65554 data
	ifstream inFile("../../data/Reaction Network (Undirected).data", ios::in);
	string lineStr;
	vector<vector<string>> strArray;
	//float array[65554][29] = { 0 };
	vector<vector<float>>array;

	if (inFile.fail())
		cout << "read failure!" << endl;

	while (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		string str;
		vector<float>data;

		getline(ss, str, ','); // leave out the first one
		while (getline(ss, str, ',')) {
			data.push_back(atof(str.c_str()));
		}
		array.push_back(data);
	}
	cout << array[0][0] << endl;
	cout << "number of instances: " << array.size() << endl;
	cout << "number of attributes: " << array[0].size() << endl;
}

void quick_sort(vector<point>& s, int dimension, int l, int r) {
	/*
	* func: 对同一维度的数据进行排序，正确性经过验证
	*/
	if (l < r) {
		int i = l, j = r;
		point x = s[l];
		while (i < j) {
			while (i < j && s[j].data[dimension] >= x.data[dimension])
				j--;

			if (i < j)
				s[i++] = s[j];

			while (i < j && s[i].data[dimension] < x.data[dimension])
				i++;

			if (i < j)
				s[j--] = s[i];

		}
		s[i] = x;
		quick_sort(s, dimension, l, i - 1);
		quick_sort(s, dimension, i + 1, r);
	}
}

void random(vector<point>& point_list, int n, int m) {
	// number of data, number of dimension
	int l = 0, r = 1000;        // number range
	srand(time(NULL));

	for (int i = 0; i < n; i++) {

		vector<int>curr_data;
		for (int j = 0; j < m; j++)
			curr_data.push_back(rand() % (r - l + 1) + l);
		point_list.push_back({ i, curr_data });
	}
}

vector<vector<int>> sort_data(vector<point>point_list){
	/**
	 * func: sort data of every dimension
	 */
	int dimension = point_list[0].data.size();
	vector<vector<int>>sort_list;
	for(int i = 0; i < dimension; ++i){
		vector<int>res(point_list.size());

		quick_sort(point_list, i, 0, point_list.size() - 1);

		for(int j = 0; j < point_list.size(); ++j)
			res[j] = point_list[j].index;

		sort_list.push_back(res);
	}
	return sort_list;
}
//! tools for comparison
Comparator* generate_comparator(){

	CircuitType type = UNI;         // 比较电路的类型
    bool verbose = true ;           // 是否输出错误信息
    unsigned long d = 2;            // field extension degree 域扩展次数

    unsigned long expansion_len = 3;    //!! maximal number of digits in a number
    SecKey* secret_key = new SecKey(context);         // Create a secret key associated with the context
    secret_key->GenSecKey();             // 生成加密所需私钥

    // 生成密钥转换矩阵
    if (expansion_len > 1){
        if (context.getZMStar().numOfGens() == 1){
            std::set<long> automVals;
            long e = 1;
            long ord = context.getZMStar().OrderOf(0);
            bool native = context.getZMStar().SameOrd(0);

            if(!native)
                automVals.insert(context.getZMStar().genToPow(0, -ord));
                
            while (e < expansion_len){
                long atm = context.getZMStar().genToPow(0, ord-e);
                //cout << "Automorphism " << -e << " is " << atm << endl;
                automVals.insert(atm);
                e <<=1;
            }
            addTheseMatrices(*secret_key, automVals);
        }
        else
            addSome1DMatrices(*secret_key);
    }

    if (d > 1)
        addFrbMatrices(*secret_key);     //might be useful only when d > 1，不知道有啥用暂且先放在这里
    
    //! 创建比较器，然后发送给Cloud1，将公钥发送给User
	
    Comparator* comparator = new Comparator(context, type, d, expansion_len, *secret_key, verbose);	

	return comparator;	
}

void divide_data(vector<point>point_list, int data_range, vector<point>&c1, vector<point>&c2){
	/**
	 * func: 将数据随机拆分成两份
	 */
	c1 = point_list;
	c2 = point_list;
	srand(time(NULL));

	for(int i = 0; i < point_list.size(); ++i){
		for(int j = 0; j < point_list[0].data.size(); ++j){
			int num = rand() % point_list[i].data[j]; //! FIXME 这里没有模数据范围，而是模原始数据，保证秘密共享拆分结果都是整数，不知道会不会有问题
			c1[i].data[j] = num;
			c2[i].data[j] -= num;
		}
	}

}

void generate_beaver_set(int n, int data_range, vector<vector<int>>&c1_list, vector<vector<int>>&c2_list){
	/**
	 * func: 生成乘法所需beaver三元组
	 * note: c1_list, c2_list都是在云中生成好了的
	 * *check: 正确性经过验证
	 */
	srand(time(NULL));
	for(int i = 0; i < n; ++i){
		int a = rand() % data_range;
		int b = rand() % data_range;
		int c = a * b;

		int c1_a = rand() % a;
		int c2_a = a - c1_a;

		int c1_b = rand() % b;
		int c2_b = b - c1_b;

		int c1_c = rand() % c;
		int c2_c = c - c1_c;

		c1_list[i] = vector<int>{c1_a, c1_b, c1_c};
		c2_list[i] = vector<int>{c2_a, c2_b, c2_c};
	}
}
void secret_share_N(vector<int>N, vector<int>&N1, vector<int>&N2){
	/**
	 * func: 将由01构成的N数组，划分为由01构成的N1，N2两个点集
	 */
	srand(time(NULL));
	for(int i = 0; i < N.size(); ++i){
		N1[i] = rand() % 1;
		N2[i] = N1[i] ^ N[i];
	}
}



