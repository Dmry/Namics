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
				vector<string> tokens;
				vector<string> line_tokens;

                line = discard_whitespace(line);

				tokens.push_back( to_string(line_no) );
				split(line,":", line_tokens);

                std::copy ( line_tokens.begin(), line_tokens.end(), std::back_inserter(tokens));

                if ( tokens.size() != 5) {
                    cerr << " Line number " << line_no << " does not contain 4 items" << endl;
                    throw ERROR::INCORRECT_TOKEN_COUNT;
                }

                m_settings.push_back(tokens);

            } else if (is_not_blank(line) and is_not_commented(line) and line == "start")
				m_settings.push_back({to_string(line_no), "start"});

        } //finished reading lines

		assert_last_token_start(line_no);

        m_file.close();

        return m_settings;
}

bool Input_parser::assert_last_token_start(size_t line_no) {
	if (m_settings.back()[1] != "start") {
		m_settings.push_back(tokenized_line{ to_string(line_no), "start", " ", " ", " "}); return false;
	} else {
		return true;
	}
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
    //return !regex_match( line, regex(R"(^\/{2,} ?)") );
    return line.substr(0, 2) != "//";
}

void Input_parser::split(const string& line, const string& delimiter, tokenized_line& target) {
    tokenized_line tokens;
    std::regex sregex("([^\\" + delimiter + "]+)");
    regex_token_iterator<string::const_iterator> end; // default end iterator
    regex_token_iterator<string::const_iterator> a(line.begin(), line.end(), sregex);
    while (a!=end)
        tokens.push_back(*a++);
    target = tokens;
}

void Input::split(const string& line, const char& delimiter, tokenized_line& target) {
	const string compatibility_delimiter = string( 1, delimiter );
	dynamic_cast<Input_parser*>(m_parser)->split(line, compatibility_delimiter, target);
}

Input::Input(string filename)
: m_parser{new Input_parser(filename)}, m_file_content{m_parser->parse()}, name{filename}
{
	KEYS.push_back("start");
	KEYS.push_back("sys");
	KEYS.push_back("mol");
	KEYS.push_back("mon");
	KEYS.push_back("alias");
 	KEYS.push_back("lat");
	KEYS.push_back("newton");
	KEYS.push_back("mesodyn");
	KEYS.push_back("cleng");
	KEYS.push_back("teng");
	KEYS.push_back("output");
	KEYS.push_back("var");
	KEYS.push_back(OutputInfo::IN_CLASS_NAME);
	KEYS.push_back("state");
	KEYS.push_back("reaction");

	parseOutputInfo();
	CheckInput();
}

std::sregex_iterator Input::NestedGroup(string& input, string& delimiter_open, string& delimiter_close) {
	std::regex rx = std::regex("(\\" + delimiter_open + ")(.*)(\\" + delimiter_close + ")");
    return sregex_iterator(input.begin(), input.end(), rx);
}

bool Input::EvenSign(string input, string delimiter_open, string delimiter_close, vector<int>& open, vector<int>& close) {
    std::sregex_iterator matches = NestedGroup(input, delimiter_open, delimiter_close);
	// position 0 is full match, 1 is first bracket, 2 is what is between brackets, 3 is closing bracket
	if (!matches->empty()) {
    	open.push_back(matches->position(1));
    	close.push_back(matches->position(3));
		EvenSign((*matches)[1].str(), delimiter_open, delimiter_close, open, close);
	}
    return !regex_match((*matches)[2].str(), regex("([\\" + delimiter_open + "\\" + delimiter_close + "])"));
}

bool Input::EvenSquareBrackets(string exp, vector<int> &open, vector<int> &close) {
    return EvenSign(exp, "[", "]", open, close);
}

bool Input::EvenBrackets(string exp, vector<int> &open, vector<int> &close) {
    return EvenSign(exp, "(", ")", open, close);
}


int Input::GetNumStarts() {
	size_t count {0};
	for (const auto& line : m_file_content)
		if (line[TOKENS::KEY] == "start") ++count;
	return count;
}

bool Input::CheckParameters(string key, string brand, int start, std::vector<std::string> &KEYS, std::vector<std::string> &PARAMETERS,std::vector<std::string> &VALUES) {
	bool success=true;
	int n_start=0;

	for (const auto& line : m_file_content) {
		if ( line[TOKENS::KEY] == string("start") ) {
			n_start++;
		} else if (n_start < start) {
			if ( line[TOKENS::KEY] == key && line[TOKENS::BRAND] == brand) {
				string parameter = line[TOKENS::PARAMETER];

				if (std::find(KEYS.begin(), KEYS.end(), parameter) != KEYS.end()) {
					if (std::count(PARAMETERS.begin(), PARAMETERS.end(), parameter) == 0) {
						parameter_value_pair[parameter] = line[TOKENS::VALUE]; success = true;
					} else {
						success=false; cout <<"After 'start' " << n_start << ", in line " << line[TOKENS::LINE] << " " << key << " property '" << parameter << "' is already defined. "<< endl;
					}

				} else  {
					success=false; cout <<"In line " << line[TOKENS::LINE] << " "  << key << " property '" << parameter << "' is unknown. Select from: "<< endl;
					for (const auto& key : KEYS) std::cout << key << endl;
				}	 
			}
		} else break;
	}

	PARAMETERS.clear(); VALUES.clear();
	for(auto& parameter_first_value_second : parameter_value_pair) {
		PARAMETERS.push_back(parameter_first_value_second.first);
		VALUES.push_back(parameter_first_value_second.second);
	}

	return success;
}

bool Input::ValidateKeys() {
	for (const auto& line : m_file_content)
		if (std::find(KEYS.begin(), KEYS.end(), line[TOKENS::KEY]) == KEYS.end()) {
			cout << line[TOKENS::KEY] << " is not valid keyword in line " << line[TOKENS::LINE] << endl;
			cout << "select one of the following:" << endl;
			for (auto& key : KEYS) cout << key << endl;
			return false;
		} else return true;
}

bool Input::CheckInput(void) {
	bool success=true;

	vector<string> options { "ana", "vtk", "kal", "pro", "vec", "pos"};

	for (const auto& line : m_file_content)
		if (line[TOKENS::KEY]=="output") {
			if (std::find(options.begin(), options.end(), line[TOKENS::BRAND]) == options.end()) { cout << "Value for output extension '" + line[TOKENS::BRAND] + "' not allowed." << endl; success = false; }
			if (std::find(KEYS.begin(), KEYS.end(), line[TOKENS::BRAND]) == KEYS.end()) KEYS.push_back( line[TOKENS::BRAND] );
		}

	success = ValidateKeys();
	parseOutputInfo();

	if (!output_info.isOutputExists()) {
		cout << "Cannot access output folder '" << output_info.getOutputPath() << "'" << endl;
		success = false;
	}

	return success;
}

void Input::parseOutputInfo() {
	for (auto& line : m_file_content) {
		if (line[TOKENS::KEY] != OutputInfo::IN_CLASS_NAME) {
			continue;
		}
		output_info.addProperty(line[TOKENS::BRAND], line[TOKENS::PARAMETER], line[TOKENS::VALUE]);
	}
}

bool Input::MakeLists(int start) {
	bool success=true;
	SysList.clear();
	LatList.clear();
	NewtonList.clear();
	MonList.clear();
	MolList.clear();
	OutputList.clear();
	MesodynList.clear();
	ClengList.clear();
	TengList.clear();
	VarList.clear();
	StateList.clear();
	ReactionList.clear();

	if (!TestNum(SysList,"sys",0,1,start)) {cout << "There can be no more than 1 'sys name' in the input" << endl; }
	if (SysList.size()==0) SysList.push_back("noname");
	if (!TestNum(LatList,"lat",1,1,start)) {cout << "There must be exactly one 'lat name' in the input" << endl; success=false;}
	if (!TestNum(NewtonList,"newton",0,1,start)) {cout << "There can be no more than 1 'newton name' in input" << endl; success=false;}
	if (NewtonList.size()==0) NewtonList.push_back("noname");
	if (!TestNum(MonList,"mon",2,1000,start)) {cout << "There must be at least one 'mon name' in input" << endl; success=false;}
	if (!TestNum(StateList,"state",0,1000,start)) {cout << "There can not be more than 1000 'state name's in input" << endl; success=false;}
	if (!TestNum(ReactionList,"reaction",0,1000,start)) {cout << "There can not be more than 1000 reaction name's in input" << endl; success=false;}
	TestNum(AliasList,"alias",0,1000,start);
	if (AliasList.size()==0) AliasList.push_back("noname");
	if (!TestNum(MolList,"mol",1,1000,start)) {cout << "There must be at least one 'mol name' in input" << endl; success=false;}
	if (!TestNum(OutputList,"output",1,1000,start)) {cout << "No output defined! " << endl;}
	if (!TestNum(MesodynList,"mesodyn",0,1,start)) {cout << "There can be no more than 1 'mesodyn' engine brand name in the input " << endl; success=false;}
	if (!TestNum(ClengList,"cleng",0,1,start)) {cout << "There can be no more than 1 'cleng' engine brand name in the input " << endl; success=false;}
	if (!TestNum(TengList,"teng",0,1,start)) {cout << "There can be no more than 1 'teng' engine brand name in the input " << endl; success=false;}
	if (!TestNum(VarList,"var",0,10,start))
	if (VarList.size()==0) VarList.push_back("noname");
	return success;
}

bool Input::TestNum(std::vector<std::string> &S, string c,int num_low, int num_high, int UptoStartNumber  ) {

	bool InList=false;
	int n_starts=0;
	for (auto& line : m_file_content) {
		if (line[1]=="start") n_starts++;
		if (c==line[1] && n_starts < UptoStartNumber ){
			InList=false;
			int S_length=S.size();
			for (int k=0; k<S_length; k++) if (line[2]==S[k]) InList=true;
			if (!InList) S.push_back(line[2]);
		}
	}
	int number=S.size();
	if (number>num_low-1 && number<num_high+1) {return true;}
	return false;
}

bool Input::IsDigit( string token )
{
    return regex_match( token, regex("^-?[\\d]+$") );
}

bool Input::is_bool( string token ) {
    return regex_match( token, regex("^[01]$|^true$|^false$") );
}

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
		cout << error << " value '" << target << "' is not allowed. Select from: " << endl;
		for (auto & option : options) cout << option << endl;
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