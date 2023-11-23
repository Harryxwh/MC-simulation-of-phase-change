#ifndef ParameterBundle_HPP
#define ParameterBundle_HPP

#include <iostream>
#include <vector>
#include <map>
#include "pugixml/pugixml.hpp"

class Parameter {
private:
	const std::string name;
	const std::string type;
	const std::string numerical_condition_type;
	bool val_set;
	
	int val_int;
	double val_dbl;
	std::string val_str;
	
	int bound_int;
	double bound_dbl;
	
	const bool print_option;
	
public:
	Parameter(const std::string name_spec, const std::string type_spec, bool print_option_spec, const std::string numerical_condition_type_spec = "N/A", double bound_spec = 0) : name(name_spec), type(type_spec), numerical_condition_type(numerical_condition_type_spec), val_set(false), val_int(3), val_dbl(3.141593), val_str("pi"), bound_int(bound_spec), bound_dbl(bound_spec), print_option(print_option_spec) {
		if (type != "int" && type != "double" && type != "string" ) {
			std::cerr << "unknown type" << std::endl;
			exit(1);
		}
		if (numerical_condition_type != "N/A" && numerical_condition_type != ">" && numerical_condition_type != "<" && numerical_condition_type != ">=" && numerical_condition_type != "<=") {
			std::cerr << "unknown numerical_condition_type" << std::endl;
			exit(1);
		}
	};
	
	~Parameter() {};
	
	bool check_val_int(int val) {
		bool check = true;
		if (numerical_condition_type == ">") check = (val > bound_int);
		if (numerical_condition_type == ">=") check = (val >= bound_int);
		if (numerical_condition_type == "<") check = (val < bound_int);
		if (numerical_condition_type == "<=") check = (val <= bound_int);
		return check;
	};
	
	bool check_val_dbl(double val) {
		bool check = true;
		if (numerical_condition_type == ">") check = (val > bound_dbl);
		if (numerical_condition_type == ">=") check = (val >= bound_dbl);
		if (numerical_condition_type == "<") check = (val < bound_dbl);
		if (numerical_condition_type == "<=") check = (val <= bound_dbl);
		return check;
	};

	std::string _name() const { return name; };
	std::string _type() const { return type; };
	bool _val_set() const { return val_set; };

	void set_val(int val_spec) {
		assert(type == "int");
		assert(val_set == false);
		if (! check_val_int(val_spec)) {
			std::cerr << " ... Error (" << name << " must be " << numerical_condition_type << " " << bound_int << ")" << std::endl;
			exit(1);
		}
		if (print_option) std::cout << '\n' << "# Parameter::set_val> Set " << name << " = " << val_spec << " as an integer" << '\n' << std::endl;
		val_int = val_spec;
		val_set = true;
	};
	
	void set_val(double val_spec) {
		assert(type == "double");
		assert(val_set == false);
		if (! check_val_dbl(val_spec)) {
			std::cerr << " ... Error (" << name << " must be " << numerical_condition_type << " " << bound_dbl << ")" << std::endl;
			exit(1);
		}
		if (print_option) std::cout << '\n' << "# Parameter::set_val> Set " << name << " = " << val_spec << " as a double prec number" << '\n' << std::endl;
		val_dbl = val_spec;
		val_set = true;
	};
	
	void set_val(std::string val_spec) {
		assert(type == "string");
		assert(val_set == false);
		if (print_option) std::cout << '\n' << "# Parameter::set_val> Set " << name << " = " << val_spec << " as a string" << '\n' << std::endl;
		val_str = val_spec;
		val_set = true;
	};
	
	int _val_int() const {
		assert(type == "int");
		return val_int;
	};
	
	double _val_dbl() const {
		assert(type == "double");
		return val_dbl;
	};
	
	std::string _val_str() const {
		assert(type == "string");
		return val_str;
	};
};

class ParameterBundle {
protected:
	std::vector<Parameter> parameters;
	std::map<std::string, Parameter*> parameters_guide;
	const bool print_option;
	
public:
	ParameterBundle(std::string input_file, bool print_option_spec = false) : print_option(print_option_spec) { init(input_file); };
		
	ParameterBundle(bool print_option_spec) : print_option(print_option_spec) {};
	
	void init(std::string input_file) {
		add_parameters();
		setup_guidemap();
		read_xml(input_file);
		check_if_any_unset();
	};

	int _value_int(std::string parameter_name) {
		std::map<std::string, Parameter*>::iterator database_entry = parameters_guide.find(parameter_name);
		if (database_entry == parameters_guide.end()) {
			std::cerr << parameter_name << " is not found" << std::endl;
			exit(1);
		} else return (database_entry->second)->_val_int();
	};
	
	double _value_dbl(std::string parameter_name) {
		std::map<std::string, Parameter*>::iterator database_entry = parameters_guide.find(parameter_name);
		if (database_entry == parameters_guide.end()) {
			std::cerr << parameter_name << " is not found" << std::endl;
			exit(1);
		} else return (database_entry->second)->_val_dbl();
	};
	
protected:
	virtual void add_parameters() {
		parameters.push_back(Parameter("ID", "string", print_option));
	};

	virtual void read_xml(std::string input_file) {
		pugi::xml_document doc;
		if ( ! doc.load_file(input_file.c_str()) ) {
			std::cerr << input_file << " is not found" << std::endl;
			exit(1);
		}
		read_xml_category(doc, "ID");
	};
	
	void read_xml_category(pugi::xml_document& doc, std::string categorry_str) {
		pugi::xml_node categoray = doc.child("ParameterData").child(categorry_str.c_str());
		for (pugi::xml_node spec: categoray.children("Parameter")) {
			if (print_option) std::cout << "# ParameterBundle::read_xml_category(...)> Parameter of the categoray \"" << categorry_str << "\" : ";
			for (pugi::xml_attribute attr: spec.attributes()) {
				set_parameter(attr);
				break;
			}
		}
	};
	
private:
	void setup_guidemap() {
		for (auto& each : parameters) parameters_guide.insert( make_pair(each._name(), &each) );
	};
	
	bool set_parameter(pugi::xml_attribute& attr) {
		bool success = false;
		pugi::xml_attribute attr_val  = attr.next_attribute();
		if (print_option) std::cout << attr.value() << " = " << attr_val.value();

		for ( auto& each : parameters ) {
			const std::string parameter_name = each._name();
			std::map<std::string, Parameter*>::iterator database_entry = parameters_guide.find(parameter_name);
			if (attr.name() == std::string("Name") && attr.value() == parameter_name) {
				success = true;
				if ( (database_entry->second)->_type() == std::string("int") ) {
					(database_entry->second)->set_val(attr_val.as_int());
				} else if ( (database_entry->second)->_type() == std::string("double") ) {
					(database_entry->second)->set_val(attr_val.as_double());
				} else {
					(database_entry->second)->set_val(attr_val.as_string());
				}
			}
		}
		if (! success) {
			std::cerr << "\"" << attr.value() << "\" is defined in neither ParameterBundle nor derived classes (error!)" << '\n' << std::endl;
			exit(0);
		}
		return success;
	};
	
	void check_if_any_unset() {
		std::vector<Parameter>::iterator itr = parameters.begin();
		while (itr != parameters.end()) {
			if (! (*itr)._val_set()) break;
			itr++;
		}
		if (itr != parameters.end()) {
			std::cerr << '\n' << "ParameterBundle::check_if_any_unset> Error; \"" << itr->_name() << "\" was not found in the input file" << std::endl;
			exit(1);
		}
	};
};
//cout
class CubicLatticeParameterBundle : public ParameterBundle {
public:
	CubicLatticeParameterBundle(std::string input_file, bool print_option_spec = false) : ParameterBundle(print_option_spec) {
		init(input_file);
	};
	
	void read_xml(std::string input_file) {
		pugi::xml_document doc;
		if ( ! doc.load_file(input_file.c_str()) ) {
			std::cerr << input_file << " is not found" << std::endl;
			exit(1);
		}
		read_xml_category(doc, "ID");
		read_xml_category(doc, "System");
		read_xml_category(doc, "Model");
	};

	void add_parameters() {
		ParameterBundle::add_parameters();
		parameters.push_back(Parameter("L0", "int", print_option, ">=", 2));
		parameters.push_back(Parameter("L1", "int", print_option, ">=", 2));
		parameters.push_back(Parameter("L2", "int", print_option, ">=", 2));
		parameters.push_back(Parameter("Seed", "int", print_option, ">=", 10000));
		parameters.push_back(Parameter("N_TemperaturePoints", "int", print_option, ">=", 2));
		parameters.push_back(Parameter("Tmax", "double", print_option, ">", 0));
		parameters.push_back(Parameter("Tmin", "double", print_option, ">", 0));
		parameters.push_back(Parameter("MaxSweepStep", "int", print_option, ">", 0));
		parameters.push_back(Parameter("n_bins", "int", print_option, ">", 1));
		parameters.push_back(Parameter("mcs_thermalization", "int", print_option, ">", 0));
		parameters.push_back(Parameter("n_samples_per_bin", "int", print_option, ">", 0));
		parameters.push_back(Parameter("mcs_interval_btwn_bins", "int", print_option, ">", 0));
	};
	
	int _L0() { return _value_int("L0"); };
	int _L1() { return _value_int("L1"); };
	int _L2() { return _value_int("L2"); };
	int _N_TemperaturePoints() { return _value_int("N_TemperaturePoints"); };
	int _Seed() { return _value_int("Seed"); };
	double _Tmin() { return _value_dbl("Tmin"); };
	double _Tmax() { return _value_dbl("Tmax"); };
	
	int _MaxSweepStep() { return _value_int("MaxSweepStep"); }

	int _n_bins() { return _value_int("n_bins"); };
	int _mcs_thermalization() { return _value_int("mcs_thermalization"); };
	int _n_samples_per_bin() { return _value_int("n_samples_per_bin"); };
	int _mcs_interval_btwn_bins() { return _value_int("mcs_interval_btwn_bins"); };
};

#endif /* ParameterBundle_HPP */
