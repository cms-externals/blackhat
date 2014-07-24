/*
 * color_algebra.h
 *
 */
#include<iostream>
#include "BH_utilities.h"
#include "constants.h"

namespace BH {

class color_index;
class color_string;
class color_constant;
class single_color_tensor;
class color_tensor;


std::ostream& operator<<(std::ostream&, const color_index& );
std::ostream& operator<<(std::ostream&, const color_string& );
std::ostream& operator<<(std::ostream&, const color_constant& );
std::ostream& operator<<(std::ostream&, const single_color_tensor& );
std::ostream& operator<<(std::ostream&, const color_tensor& );


class color_index {
	protected:
		enum index_type {adjoint,fundamental,unknown};
		index_type kind;
		size_t label;
		friend std::ostream& operator<<(std::ostream&, const color_index& );
	public:
		// constructors
		color_index();
		color_index(std::string ct,int lb);
		friend class color_string;
		friend class single_color_tensor;
		friend class color_tensor;
		// destructor
//		~color_index();
};

class adjoint_index: public color_index {
	public:
		// constructor
		adjoint_index(int lb);
};

class fundamental_index: public color_index {
	public:
		// constructor
		fundamental_index(int lb);
};


class color_string {
	private:
		// can be any natural number
		size_t adjoint_degree;
		std::vector<size_t> adjoint_color_labels;
		std::vector<color_index> adjoint_color_indices;
		// can be either 0 or 1
		size_t fundamental_degree;
		std::vector<size_t> fundamental_color_labels;
		std::vector<color_index> fundamental_color_indices;
		friend	bool operator==(const color_string& , const color_string&);
		friend	bool operator<(const color_string& , const color_string&);
		friend	bool compare_cs(color_string* , color_string*);
		friend std::ostream& operator<<(std::ostream&, const color_string& );
	public:
		// constructors
		color_string();
		// if fundamental indices appear, should be in the beginning and end of the vector! --- checked in the constructor
		color_string(std::vector<color_index*> cis);
		color_string(const color_index& ci1,const color_index& ci2,
				const color_index& ci3=color_index(),const color_index& ci4=color_index(),const color_index& ci5=color_index(),
				const color_index& ci6=color_index(),const color_index& ci7=color_index(),const color_index& ci8=color_index());
		// copy constructor
		color_string(const color_string& scs);
		//return minus sign when odd number of adjoint matrices	
		int conjugate();
		//simplification functions
		color_tensor* simplify(const color_string&);
		color_tensor* simplify();
//		std::vector<color_constant > join_prefactors(const color_string&, const color_string);
		color_string& operator=(const color_string& scs);
		bool rotate_trace();
		friend class single_color_tensor;
		friend class color_tensor;
};

bool operator==(const color_string& cs1, const color_string& cs2);
bool operator<(const color_string& cs1, const color_string& cs2);
bool compare_cs(color_string* cs1, color_string* cs2);
	

// it reprensents a polinomial on a given parameter --- it can be over several bases
// in principle it represents Nc polynomials --- it might be adapted in the future for combinations like (Nc^2-1)
class color_constant {
	std::vector<multi_precision_constant> positive_coeffs;
	std::vector<multi_precision_fraction> positive_coeffs_fracs;
	std::vector<multi_precision_constant> negative_coeffs;
	std::vector<multi_precision_fraction> negative_coeffs_fracs;
	// color_constant operations
	friend const color_constant operator+(const color_constant& cc1,const color_constant& cc2);
	friend const color_constant operator-(const color_constant& cc1,const color_constant& cc2);
	friend const color_constant operator*(const color_constant& cc1,const color_constant& cc2);
	friend const color_constant operator*(int i,const color_constant& cc2);
	friend bool operator==(const color_constant& cc1, const color_constant& cc2);
	friend std::ostream& operator<<(std::ostream&, const color_constant& );
	public:
		// constructors
		color_constant(){};
		color_constant(int num,int power);
		color_constant(int num,int den,int power);
		bool is_zero();
		double eval();
		void project_to_Nc_powers(int,int);
		color_constant& operator+=(const color_constant& cc1);
		color_constant& operator=(const color_constant& cc1);
		friend class single_color_tensor;
		friend class color_tensor;
};

bool operator==(const color_constant& cc1, const color_constant& cc2);

// the (external) product of single_color_tensor is a single_color_tensor
// the sum of two single_color_tensor is not always a single_color_tensor, but a color_tensor in general
class single_color_tensor {
	private:
	color_constant color_coeff;
	std::vector<color_string* > color_string_vector;
	friend const single_color_tensor operator*(const single_color_tensor& cc1,const single_color_tensor& cc2);
	friend const single_color_tensor operator*(const color_constant& cc,const single_color_tensor& cc2);
	friend bool operator==(const single_color_tensor&, const single_color_tensor&);
	friend bool operator<(const single_color_tensor&, const single_color_tensor&);
	friend bool compare_sct(single_color_tensor*, single_color_tensor*);
	friend std::ostream& operator<<(std::ostream&, const single_color_tensor& );
	public:
	// constructors
		single_color_tensor(){};
		//copy constructor
		single_color_tensor(const single_color_tensor& sct);
		single_color_tensor(const color_constant& coeff, const color_string& csv);
		single_color_tensor(const color_constant& coeff, std::vector<color_string* > csv);
		single_color_tensor(const color_constant& coeff, 
				color_index ci1,
				color_index ci2,
				color_index ci3=color_index(),
				color_index ci4=color_index(),
				color_index ci5=color_index(),
				color_index ci6=color_index(),
				color_index ci7=color_index(),
				color_index ci8=color_index(),
				color_index ci9=color_index(),
				color_index ci10=color_index(),
				color_index ci11=color_index());
		color_tensor* simplify();
		bool is_zero();
		bool sort();
		void conjugate();
		bool sortedQ();
		color_constant* project_to_color_constant();
		color_constant* project_to_single_color_tensor(const single_color_tensor& sct );

		single_color_tensor& operator=(const single_color_tensor& sct);
		friend class color_tensor;
};


bool operator==(const single_color_tensor& sct, const single_color_tensor& sct2);
bool operator<(const single_color_tensor& sct, const single_color_tensor& sct2);
bool compare_sct(single_color_tensor* sct, single_color_tensor* sct2);



class color_tensor {
	private:
	std::vector<single_color_tensor* > single_color_tensors;
	// color_tensor operations
	friend const color_tensor operator+(const color_tensor& cc1,const color_tensor& cc2);
	friend const color_tensor operator-(const color_tensor& cc1,const color_tensor& cc2);
	friend const color_tensor operator*(const color_tensor& cc1,const color_tensor& cc2);
	friend const color_tensor operator*(const color_constant& cc,const color_tensor& cc1);
	friend std::ostream& operator<<(std::ostream&, const color_tensor& );
	public:
		// constructors
	color_tensor(){};
	color_tensor(const single_color_tensor& sct);
	color_tensor(std::vector<single_color_tensor*>  v_sct);
	color_tensor(const color_constant& coeff, std::vector<color_string* > csv);
	color_tensor(const color_constant& coeff, const color_string& csv);
	color_tensor(const color_constant& coeff, 
				color_index ci1,
				color_index ci2,
				color_index ci3=color_index(),
				color_index ci4=color_index(),
				color_index ci5=color_index(),
				color_index ci6=color_index(),
				color_index ci7=color_index(),
				color_index ci8=color_index(),
				color_index ci9=color_index(),
				color_index ci10=color_index(),
				color_index ci11=color_index());
	color_tensor& operator=(const color_tensor& ct);
	void conjugate();
	bool single_simplify();
	void sort();
	bool shorten();
	void simplify();
	bool is_zero();
	color_constant* project_to_color_constant();
	color_constant* project_to_color_tensor(const color_tensor& ct);
};



}
