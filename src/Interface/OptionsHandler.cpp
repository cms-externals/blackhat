/*
 * OptionsHandler.cpp
 *
 *  Created on: 26 Mar 2010
 *      Author: daniel
 */

#include "OptionsHandler.h"
#include "OptionsHandler.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

namespace BH {

namespace Tools {

singleValueOption::singleValueOption(const std::string& name,const std::string& value,const std::string& helpString): NamedOption(name,helpString), d_value(value) {};

bool singleValueOption::process(std::istream& is,std::string& message) const {
	Parser<std::string> parser;
	string answer=parser.parse(is);
	if ( answer == d_value  ){
		message = "" ;
		return true;
	} else {
		//if we get there, we didn't undersstand the answer
		message = "ERROR: The only possible value is: \'" + d_value + "\' for option " + getName() + "." ;
		return false;
	}
}

std::ostream& singleValueOption::printHelp(std::ostream& os) const {
	return os << d_helpString << std::string(" The only available value is: ") << d_value << ".";
}


bool yesOrNoOption::process(std::istream& is,std::string& message) const {
	string answer;
	is >> answer;
	if ( answer == "Yes" || answer == "On" || answer == "yes" || answer == "on" || answer == "YES" || answer == "ON"  ){
		d_valueToSet=true;
		message = "" ;
		return true;
	}
	if ( answer == "No" || answer == "Off" || answer == "no" || answer == "off" || answer == "NO" || answer == "OFF"  ){
		d_valueToSet=false;
		message = "" ;
		return true;
	}
	//if we get there, we didn't undersstand the answer
	message = "ERROR: could not understand the answer: \'" + answer + "\' for option " + getName() + "." ;
	return false;
}

yesOrNoOption::yesOrNoOption(const std::string& name,bool& valueToSet,const std::string& helpString): NamedOption(name,helpString), d_valueToSet(valueToSet) {};

AlwaysErrorOption::AlwaysErrorOption(const std::string& name,const std::string& errorMessage ,const std::string& helpString): NamedOption(name,helpString), d_errorMessage(errorMessage) {};

bool AlwaysErrorOption::process(std::istream& is,std::string& message) const {
	message=d_errorMessage;
	return false;
} ;

std::ostream& AlwaysErrorOption::printHelp(std::ostream& os){
	return os << d_helpString << " NOTE: " << d_errorMessage ;
};

multipleValueOptionWithTableUpdate::multipleValueOptionWithTableUpdate(const std::string& name,
		const std::string& s1,
		const std::string& helpString,setable* ST): NamedOption(name,helpString), d_ST(ST) {
	d_names.push_back(s1);
};
multipleValueOptionWithTableUpdate::multipleValueOptionWithTableUpdate(const std::string& name,
		const std::string& s1,
		const std::string& s2,
		const std::string& helpString,setable* ST): NamedOption(name,helpString), d_ST(ST) {
	d_names.push_back(s1);
	d_names.push_back(s2);

};
multipleValueOptionWithTableUpdate::multipleValueOptionWithTableUpdate(const std::string& name,
		const std::string& s1,
		const std::string& s2,
		const std::string& s3,
		const std::string& helpString,setable* ST): NamedOption(name,helpString), d_ST(ST) {
	d_names.push_back(s1);
	d_names.push_back(s2);
	d_names.push_back(s3);

};
multipleValueOptionWithTableUpdate::multipleValueOptionWithTableUpdate(const std::string& name,
		const std::string& s1,
		const std::string& s2,
		const std::string& s3,
		const std::string& s4,
		const std::string& helpString,setable* ST): NamedOption(name,helpString), d_ST(ST) {
	d_names.push_back(s1);
	d_names.push_back(s2);
	d_names.push_back(s3);
	d_names.push_back(s4);

};
multipleValueOptionWithTableUpdate::multipleValueOptionWithTableUpdate(const std::string& name,
		const std::string& s1,
		const std::string& s2,
		const std::string& s3,
		const std::string& s4,
		const std::string& s5,
		const std::string& helpString,setable* ST): NamedOption(name,helpString), d_ST(ST) {
	d_names.push_back(s1);
	d_names.push_back(s2);
	d_names.push_back(s3);
	d_names.push_back(s4);
	d_names.push_back(s5);

};

bool multipleValueOptionWithTableUpdate::process(std::istream& is,std::string& message) const {
	std::string name;
	Parser<std::string> parser;
	name = parser.parse(is);
	std::vector<std::string>::const_iterator beg = d_names.begin();
	std::vector<std::string>::const_iterator end = d_names.end();
	std::vector<std::string>::const_iterator it = find(d_names.begin(),d_names.end(),name);

	if  ( it != end ){
		d_ST->set(getName(),name) ;
		message = "" ;
		return true;
	}
	else {
		message = "ERROR: unsupported " + getName() + " value \'" + name + "\'. The only supported values are \n   " ;
		std::vector<std::string>::const_iterator it;
		for (it = beg; it != end;++it ){
			message += (*it) ;
			message += "\n   " ;
		}
		return false;
	}
}

void OptionsHandler::printHelp(ostream& os) const {
	map<string,option*>::const_iterator it = d_optionMap.begin();

	while (it != d_optionMap.end() ){
		os << (*it).first << ": \t";
		(*it).second->printHelp(os);
		os << "\n\n";
		++it;
	}

};

bool OptionsHandler::process(std::istream& is,std::string& message) const {
	string line;
	getline(is,line);
	if ( d_debug ){
		cout << "treating option line \'" << line << "\'" << endl;
	}
	stringstream ss(line);
	string optionName ;
	Parser<std::string> parser;
	optionName = parser.parse(ss);

	const map<string,option*>::const_iterator it = d_optionMap.find(optionName);

	if (it != d_optionMap.end() ){
		bool result=(*it).second->process(ss,message);
		if (result){
			d_state=success;
			return true;
		} else {
			d_state=failed;
			return false;
		}
	} else {
		d_state=unknown;
		message= "No option called " + optionName + " found.";
		return false;
	}

}

void OptionsHandler::add(option* opt){
	d_optionMap.insert(make_pair(opt->getName(),opt));
}

OptionsHandler::~OptionsHandler(){
	std::map<std::string,option*>::iterator it=d_optionMap.begin();
	std::map<std::string,option*>::iterator end=d_optionMap.end();
	for (;it!=end;it++){
		delete (*it).second;
	}
}

template class multipleValueOption<int>;
template class multipleValueOption<double>;
template class multipleValueOption<string>;

template class ValueSettingOption<double>;
template class ValueSettingOption<int>;

template class SettingsTableOption<double>;
template class SettingsTableOption<int>;


}

}
