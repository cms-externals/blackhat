/*
 * OptionsHandler.h
 *
 *  Created on: 26 Mar 2010
 *      Author: daniel
 */

#ifndef OPTIONSHANDLER_H_
#define OPTIONSHANDLER_H_

#include <string>
#include <map>
#include <vector>
#include <iosfwd>

namespace BH {

class setable;

namespace Tools {

struct ParserError  {
	std::string d_str;
	ParserError(std::string str): d_str(str) {};
};

struct ParserTypeError : ParserError {
	ParserTypeError(const std::string& str): ParserError(str){};
};

struct ParserNoInput  {
	ParserNoInput(){};
};

template <class T> class Parser {
public :
	Parser(){};
	T parse(std::istream& is) const;
};



class option {
public:
	virtual bool process(std::istream& is,std::string& os) const =0;
	virtual std::ostream& printHelp(std::ostream&) const = 0;
	virtual const std::string& getName() const =0 ;
	virtual ~option(){};
};



class NamedOption : public option {
protected:
	std::string d_optionName;
	std::string d_helpString;
public :
	NamedOption(const std::string& name,const std::string& helpString): d_optionName(name), d_helpString(helpString){};
	virtual std::ostream& printHelp(std::ostream& os) const {
		return os << d_helpString;
	}
	virtual const std::string& getName() const {
		return  d_optionName;
	}
	virtual ~NamedOption(){}
};


class yesOrNoOption : public NamedOption {
	bool& d_valueToSet;
public :
	yesOrNoOption(const std::string& name,bool& valueToSet,const std::string& helpString);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual ~yesOrNoOption(){}
};

class AlwaysErrorOption : public NamedOption {
	std::string d_errorMessage;
public :
	AlwaysErrorOption(const std::string& name,const std::string& errorMessage ,const std::string& helpString);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual std::ostream& printHelp(std::ostream&);
	virtual ~AlwaysErrorOption(){}
};


template <class TargetType> class multipleValueOption : public NamedOption {
	std::map <std::string,TargetType> d_nameMap;
	TargetType& d_valueToSet;
public :
	multipleValueOption(const std::string& name,TargetType& valueToSet,
			const std::string&, const TargetType&,
			const std::string& helpString);
	multipleValueOption(const std::string& name,TargetType& valueToSet,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string& helpString);
	multipleValueOption(const std::string& name,TargetType& valueToSet,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string& helpString);
	multipleValueOption(const std::string& name,TargetType& valueToSet,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string& helpString);
	multipleValueOption(const std::string& name,TargetType& valueToSet,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string&, const TargetType&,
			const std::string& helpString);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual ~multipleValueOption(){}
};

class multipleValueOptionWithTableUpdate : public NamedOption  {
	std::vector<std::string> d_names;
	setable* d_ST;
public :
	multipleValueOptionWithTableUpdate(const std::string& name,
			const std::string& s1,
			const std::string& helpString,setable* ST);
	multipleValueOptionWithTableUpdate(const std::string& name,
			const std::string& s1,
			const std::string& s2,
			const std::string& helpString,setable* ST);
	multipleValueOptionWithTableUpdate(const std::string& name,
			const std::string& s1,
			const std::string& s2,
			const std::string& s3,
			const std::string& helpString,setable* ST);
	multipleValueOptionWithTableUpdate(const std::string& name,
			const std::string& s1,
			const std::string& s2,
			const std::string& s3,
			const std::string& s4,
			const std::string& helpString,setable* ST);
	multipleValueOptionWithTableUpdate(const std::string& name,
			const std::string& s1,
			const std::string& s2,
			const std::string& s3,
			const std::string& s4,
			const std::string& s5,
			const std::string& helpString,setable* ST);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual ~multipleValueOptionWithTableUpdate(){}
};


template <class TargetType> class ValueSettingOption : public NamedOption {
	TargetType& d_valueToSet;
public :
	ValueSettingOption(const std::string& name,TargetType& valueToSet,
			const std::string& helpString);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual ~ValueSettingOption(){}
};

template <class TargetType> class SettingsTableOption : public NamedOption {
	setable* d_ST ;
public :
	SettingsTableOption(const std::string& name,setable* ST,
			const std::string& helpString);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual ~SettingsTableOption(){}
};


class singleValueOption : public NamedOption {
	std::string d_value;
public :
	singleValueOption(const std::string& name,const std::string& value,
			const std::string& helpString);
	virtual bool process(std::istream& is,std::string& message) const ;
	virtual std::ostream& printHelp(std::ostream& os) const;
	virtual ~singleValueOption(){}
};


class OptionsHandler {
	bool d_debug;
	std::map<std::string,option*> d_optionMap;

public:
	enum returnCode { success=0, unknown=1, failed=2 };
	mutable returnCode d_state;
	OptionsHandler() : d_debug(false) {};
	void printHelp(std::ostream& os) const;
	bool process(std::istream& is,std::string& message) const ;
	void enableDebug(){ d_debug = true; };
	void disableDebug(){ d_debug = false; };
	void add(option* opt);
	virtual ~OptionsHandler();
};



} /* Tools */

} /* BH */


#endif /* OPTIONSHANDLER_H_ */
