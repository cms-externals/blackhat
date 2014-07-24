namespace BH {


inline const process& cutD::get_process(size_t cor) const {
	return process_list[cor-1];
}

inline long cutD::get_process_code(size_t cor) const {
	return _process_code[cor-1];
}



inline size_t part::nc() const {
	return corner.size();
}

inline size_t part::np() const {
	int res=0;
	for (size_t i=0; i< nc(); i++){
		res+=corner[i].size();
	}
	return res;
}

inline const std::vector<plabel>& part::c(int i) const {
	return corner[i-1];
}
inline void part::set_code(int newcode){
	code=newcode;
}
inline int part::get_code() const {
	return code;
}

inline size_t raw_part::nc() const {
	return _corner.size();
}

inline const std::vector<particle*>& raw_part::c(int i) const {
	return _corner[i-1];
}
inline const std::vector<particle*>& raw_part::corner(int i) const {
	return _corner[i-1];
}

inline const std::vector<size_t>& raw_part::ind(int i) const {
	return _indices[i-1];
}

inline const std::vector<size_t>& raw_part::indices(int i) const {
	return (_indices[i-1]);
}

inline void raw_part::set_code(int newcode) {
	_code=newcode;
}

inline int raw_part::get_code() const {
	return _code;
}

}
