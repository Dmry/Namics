#include "input.h"

typedef IInput_parser::tokenized_line tokenized_line;

Input_parser::Input_parser(string filename)
: m_file(filename)
{
    if (!m_file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(0);
    }
}

vector<tokenized_line> Input_parser::parse() {
	string line;
	size_t line_no = 1;

	for (;getline(m_file, line); ++line_no) {
		if ( is_not_blank(line) and is_not_commented(line) and is_not_start(line)) {
			vector<string> tokens, line_tokens;

			line = discard_whitespace(line);

			tokens.push_back( to_string(line_no) );
			split(line,':', line_tokens);

			// Copy the line as separate strings (tokens) into the tokens vector
			std::copy ( line_tokens.begin(), line_tokens.end(), std::back_inserter(tokens));

			// Verify we have a line number and 4 tokens
			if ( tokens.size() != 5) {
				cerr << " Line number " << line_no << " does not contain 4 items" << endl;
				throw ERROR::INCORRECT_TOKEN_COUNT;
			}

			// Add the tokenized line to the end of the vector that will be returned
			m_settings.push_back(tokens);


		// If the line is commented or blank, do nothing, but if it's start: add start plus line number to the returned vector.
		} else if (is_not_blank(line) and is_not_commented(line) and line == "start")
			m_settings.push_back({to_string(line_no), "start", " ", " ", " "});

	} //finished reading lines

	assert_last_token_start(line_no);

	m_file.close();

	return m_settings;
}

void Input_parser::assert_last_token_start(size_t line_no) {
	if (m_settings.back()[1] != "start")
		m_settings.push_back(tokenized_line{ to_string(line_no), "start", " ", " ", " "});
}

bool Input_parser::is_not_start(const string& line) {
        string START = "start";
        return line != START;
}

string Input_parser::discard_whitespace(const string& line) {
    regex whitespace("[\\s]+");
    return regex_replace(line,whitespace,"");
}

bool Input_parser::is_not_blank(const string& token) {
    return !regex_match( token, regex(R"(^[\\s]{0,}$)") );
}

bool Input_parser::is_not_commented(const string& line) {
    return line.substr(0, 2) != "//";
}

void Input_parser::split(const string& line, const char& delimiter, tokenized_line& target) {
    std::istringstream iss { line };
	string token;
  	while ( std::getline( iss, token, delimiter ) )
	  	target.push_back(token);
}

void Input::split(const string& line, const char& delimiter, tokenized_line& target) {
	dynamic_cast<Input_parser*>(m_parser)->split(line, delimiter, target);
}

Input::Input(string filename)
: m_parser{new Input_parser(filename)}, m_file_content{m_parser->parse()}, name{filename},
out_options{ "ana", "vtk", "kal", "pro", "vec", "pos"}
{
	// Listmap: A map that can be accessed using the name of the module.
	// Rangemap: A map that tells the input verification how many of these guys to expect.
	ListMap["sys"] = &SysList;			RangeMap["sys"] = Require_input_range(0,1);
	ListMap["mol"] = &MolList;			RangeMap["mol"] = Require_input_range(1,1000);
	ListMap["mon"] = &MonList;			RangeMap["mon"] = Require_input_range(2,1000);
	ListMap["alias"] = &AliasList;		RangeMap["alias"] = Require_input_range(0,1000);
 	ListMap["lat"] = &LatList;			RangeMap["lat"] = Require_input_range(1,1);
	ListMap["newton"] = &NewtonList;	RangeMap["newton"] = Require_input_range(0,1);
	ListMap["mesodyn"] = &MesodynList;	RangeMap["mesodyn"] = Require_input_range(0,1);
	ListMap["cleng"] = &ClengList;		RangeMap["cleng"] = Require_input_range(0,1);
	ListMap["teng"] = &TengList;		RangeMap["teng"] = Require_input_range(0,1);
	ListMap["output"] = &OutputList;	RangeMap["output"] = Require_input_range(1,1000);
	ListMap["var"] = &VarList;			RangeMap["var"] = Require_input_range(0,10);
	ListMap["state"] = &StateList;		RangeMap["state"] = Require_input_range(0,1000);
	ListMap["reaction"] = &ReactionList;RangeMap["reaction"] = Require_input_range(0,1000);

	// We also want to recognize starts though.
	KEYS.push_back("start");

	// Add all the keys of the maps above to the KEYS vector for compatibility.
	for (auto& map : ListMap)
		KEYS.push_back(map.first);

	for (auto& option : out_options)
		KEYS.push_back(option);

	// So we can use binary search!
	std::sort(KEYS.begin(), KEYS.end());

	PreProcess();
	parseOutputInfo();
	CheckInput();
}

bool Input::EvenDelimiters(string exp, vector<int>& open, vector<int>& close, char delimiter_open, char delimiter_close) {
	std::stack<char> S;
	for (size_t i = 0; i < exp.size(); i++) {
		if (exp[i] == delimiter_open ) {
			S.push(exp[i]);
			open.push_back(i);
		}
		else if (exp[i] == delimiter_close) {
			close.push_back(i);
			if (S.empty() or exp[i] != delimiter_close)
				return false;
			else
				S.pop();
		}
	}
	return S.empty();
}

bool Input::EvenSquareBrackets(string exp,vector<int> &open, vector<int> &close) {
	return EvenDelimiters(exp, open, close, '[', ']');
}

bool Input::EvenBrackets(string exp,vector<int> &open, vector<int> &close) {
	return EvenDelimiters(exp, open, close, '(', ')');
}

int Input::GetNumStarts() {
	return m_settings.size();
}

void Input::PreProcess() {
	int start=1;
	size_t errors{0};

	for (const auto line : m_file_content) {
		if ( line[KEY] == string("start") and line[LINE] != m_file_content.back()[LINE] ) {
			m_settings[start+1] = m_settings[start];
			start++;
		} else {
			if (m_settings[start][line[KEY]][line[BRAND]].find(line[PARAMETER]) != m_settings[start][line[KEY]][line[BRAND]].end()
				and find(out_options.begin(), out_options.end(), line[PARAMETER]) != out_options.end()) {
					cerr << LineError(line) << "is already defined. "<< endl;
					++errors;
			} else {
					m_settings[start][line[KEY]][line[BRAND]][line[PARAMETER]] = line[VALUE];
			}
		}
	}

	if (errors > 0) {
		std::cerr << errors << " errors in input, please correct!" << endl;
		throw INPUT_ERROR::ALREADY_DEFINED;
	}

	ValidateBrands();
	ValidateInputKeys();
}

void Input::ValidateBrands() {
	// Validates number of brands in the input file and checks if they're consistent with the requested output
	for (int start = 1 ; start <= GetNumStarts() ; ++start) {
		if (not MakeLists(start))
			throw INPUT_ERROR::INPUT_RANGE;

		for (auto& output_option : out_options)
			if (OutputTypeRequested(start,output_option))
				for (auto& key : m_settings[start][output_option]) {
					AssertRequestedOutputKeyExists(start, key.first);
					for (auto& brand : key.second)
						AssertRequestedOutputBrandExists(start, key.first, brand.first);
				}
	}
}

void Input::ValidateInputKeys() {
	for (auto& keys : m_settings[GetNumStarts()])
		if (not ValidateKey(KEYS, keys.first)) {
			cerr << "Unknown key '" << keys.first << "' in input file!" << endl;
			throw INPUT_ERROR::UNKNOWN;
		}
}

bool Input::OutputTypeRequested(int start, const string& output_option) {
	return m_settings[start].find(output_option) !=  m_settings[start].end();
}

void Input::AssertRequestedOutputKeyExists(int start, const std::string& requested_key) {
	if (not std::binary_search(KEYS.begin(), KEYS.end(), requested_key)) {
		cerr << "In start " << start << ": Requested output for unrecognized module '" << requested_key << "'. Please select from:" << endl << PrintList(KEYS);
		throw INPUT_ERROR::UNKNOWN;
	}
}

void Input::AssertRequestedOutputBrandExists(int start, const std::string& requested_key, const string& requested_brand) {
	auto& list = *ListMap[requested_key];
	if (requested_brand != "*")
		// Find if brand for requested output exists or not
		if (std::find(list.begin(), list.end(), requested_brand) == list.end()) {
			std::cerr << "In start " << start << ": Name '" << requested_brand << "' " << "for " << requested_key << " in output settings not recognized. Please select from:" << endl << PrintList(list);
			throw INPUT_ERROR::UNKNOWN;
		}
}

string Input::LineError(const tokenized_line& line) {
	ostringstream error_msg;
	error_msg << "In line " << line[LINE] << " " << line[KEY] << " : " << line[BRAND] << " property '" << line[PARAMETER] << "' ";
	return error_msg.str();
}

string Input::PrintList(const std::vector<std::string>& list) {
	ostringstream error_msg;
	for (const auto& key : list) error_msg << key << endl;
	return error_msg.str();
}

bool Input::ValidateKey(const std::vector<std::string>& KEYS_, const std::string& key) {
	return std::binary_search(KEYS_.begin(), KEYS_.end(), key);
}

bool Input::CheckParameters(string key, string brand, int start, std::vector<std::string> &KEYS_, std::vector<std::string> &PARAMETERS,std::vector<std::string> &VALUES) {
	bool success=true;

	//For faster searching and to make sure people don't rely on the order of vectors...
	std::sort(KEYS_.begin(), KEYS_.end());

	PARAMETERS.clear(); VALUES.clear();
	for(const auto& key_value : m_settings[start][key][brand] ) {
		if (ValidateKey(KEYS_, key_value.first)) {
			PARAMETERS.emplace_back(key_value.first);
			VALUES.emplace_back(key_value.second);
		}
		else {
			std::cerr << key << " : " << brand << " property '" << key_value.first << "' is unknown." << endl;
			success = false;
		}
	}

	if (success == false)
		std::cerr << "Please select from:" << endl << PrintList(KEYS_);

	return success;
}

bool Input::CheckInput(void) {
	bool success=true;

	for (auto& start : m_settings)
		for (const auto & output_brand : start.second["output"]   ) {
			if (std::find(out_options.begin(), out_options.end(), output_brand.first) == out_options.end()) {
				cout << "Value for output extension '" + output_brand.first + "' not allowed." << endl;
				success = false;
			}
		}

	parseOutputInfo();

	if (!output_info.isOutputExists()) {
		cout << "Cannot access output folder '" << output_info.getOutputPath() << "'" << endl;
		throw INPUT_ERROR::OUTPUT_FOLDER;
	}

	return success;
}

void Input::parseOutputInfo() {
	for (auto& line : m_file_content) {
		if (line[KEY] != OutputInfo::IN_CLASS_NAME) {
			continue;
		}
		output_info.addProperty(line[BRAND], line[PARAMETER], line[VALUE]);
	}
}

bool Input::MakeLists(int start) {
	bool result = true;

	for (auto& map : ListMap) {
		//map: second: List vector, first: String name
		//rangemap: second : high value, first: low value
		map.second->clear();
		if (!LoadList( *map.second, map.first, RangeMap[map.first].first, RangeMap[map.first].second, start))
			result = false;
	}

	ForceNonEmptyNameFor("sys", "newton", "alias", "var");

	if (OutputList.empty()) {
		std::cout << "WARNING: No output defined for start " << start << "! " << std::endl;
	}

	return result;
}

void Input::ForceNonEmptyNameFor(const std::string name) {
	auto it = ListMap.find(name);
	if (it != ListMap.end() and it->second->empty())
		it->second->push_back("noname");
}

void Input::ForceNonEmptyNameFor(const std::string first_name, const std::string other_names...) {
	auto it = ListMap.find(first_name);
	if (it != ListMap.end() and it->second->empty())
		it->second->push_back("noname");
	ForceNonEmptyNameFor(other_names);
}

bool Input::LoadList(std::vector<std::string> &list, string key, int num_low, int num_high, int start  ) {

	for (auto& brands : m_settings[start][key])
		list.push_back(brands.first);
	
	return VerifyRange(list.size(), num_low, num_high, key);
}

bool Input::VerifyRange(int number, int num_low, int num_high, std::string& key) {
	if (number >= num_low and number <= num_high) {
		return true;
	} else {
		if (num_low != num_high)
			std::cerr << "Only between " << num_low << " and " << num_high << " '" << key << "' names are allowed in the input file." << endl;
		else
			if (num_low == 1)
				std::cerr << "There must be exactly one '" << key << "' name in the input file." << endl;
			else 
				std::cerr << "Only " << num_low << " '" << key << "' names are allowed in the input file." << endl;
		return false;
	}
}

bool Input::IsDigit( string token )
{
    return regex_match( token, regex("^-?[\\d]+$") );
}

bool Input::is_bool( string token ) {
    return regex_match( token, regex("^[01]$|^true$|^false$") );
}


/****** All of the below is added for compatibility with older code, but I advise using the to_value template from the header ******/

int Input::Get_int(string input, int default_value) {
	try { return stoi(input); }
	catch (logic_error& error ) { return default_value; }
}

bool Input::Get_int(string input, int &target, const std::string &error) {
	try { target = stoi(input); return true; }
	catch (logic_error& stoi_error) { cout << error << endl; return false; }
}

bool Input::Get_int(string input, int &target, int low, int high, const std::string &error) {
	Get_int(input, target, error);
	if (target<low || target>high) {cout << "Value out of range: " << error << endl ; return false; }
	else return true;
}

string Input::Get_string(string input, const string &target) {
	if (!input.empty()) { return input;} else {return target;}

}
bool Input::Get_string(string input, string &target, const  std::string &error) {
	if (!input.empty()) { target=input; return true;}
	else { cout << error << endl; return false;}
}
bool Input::Get_string(string input, string &target, std::vector<std::string>&options, const std::string &error) {
	if (!input.empty()) target=input;
	else {cout << error << endl; return false;}

	if (!InSet(options,target)) {
		cout << error << " value '" << target << "' is not allowed. Select from: " << endl << PrintList(options);
		return false;
	} else return true;
}

Real Input::Get_Real(string input, Real default_value) {
	try { return stod(input); }
	catch (logic_error& error ) { return default_value; }
}
bool Input::Get_Real(string input, Real &target, const std::string &error) {
	try { target = stod(input); return true; }
	catch (logic_error& stoi_error) { cout << error << endl; return false; }
}
bool Input::Get_Real(string input, Real &target, Real low, Real high, const  std::string &error) {
	Get_Real(input, target, error);
	if (target < low or target > high) {cout << "Value out of range: " << error << endl ; return false; }
	else return true;
}

bool Input::Get_bool(string input, bool default_value) {
	bool output = false;
	if( is_bool(input) )
		return to_value(input, output);
	else return default_value;
}

bool Input::Get_bool(string input, bool &target, const std::string &error) {
	if( is_bool(input) ) { to_value(input, target); return true; }
	else { cout << error << endl; return false; }
}

bool Input:: LoadItems(const int start, const string key,std::vector<std::string> &Out_key, std::vector<std::string> &Out_name, std::vector<std::string> &Out_prop) {

	brand_settings_map output;
	string wildcard_char = "*";

	try {
		output = m_settings[start].at(key);
	} catch (out_of_range) {
		cerr << "No output found for '" << key << "'. This is definitely a developer's error." <<  endl;
		throw INPUT_ERROR::OUT_NAME_NOT_FOUND;
	}

	for (auto& keys : output) {
		auto& requested_key = keys.first;
		for (auto& brands : keys.second) {
			auto& requested_brand = brands.first;
			if (requested_brand == wildcard_char)
				for (auto& brand : *ListMap[requested_key]) {
					Out_key.push_back(requested_key); Out_name.push_back(brand); Out_prop.push_back(brands.second);
				}
			else
				{
					Out_key.push_back(requested_key); Out_name.push_back(requested_brand); Out_prop.push_back(brands.second);
				}
					
		}
	}
	return true;
}

// Seems to be used by Lattice
bool Input::ReadFile(string fname, string &In_buffer) {
	ifstream this_file;
	bool success=true;
	bool add;
	this_file.open(fname.c_str());
	std:: string In_line;
	if (this_file.is_open()) {
		while (this_file) {
			add=true;
			std::getline(this_file,In_line);
			In_line.erase(std::remove(In_line.begin(), In_line.end(), ' '), In_line.end());
			if (In_line.length()==0) add = false;
			if (In_line.length()>2) {if (In_line.substr(0,2) == "//") {add = false;}}
			if (add) {In_buffer.append(In_line).append("#");};
		}
		this_file.close();
		if (In_buffer.size()==0) {cout << "File " + fname + " is empty " << endl; success=false; }
	} else {cout <<  "Inputfile " << fname << " is not found. " << endl; success=false; }
	return success;
}
