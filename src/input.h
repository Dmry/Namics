#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <algorithm>
#include <vector>
#include <utility>
#include <regex>
#include <iostream>
#include <exception>
#include <map>
#include <deque>
#include <stack>
#include <functional>
#include <fstream>
#include "output_info.h"

#define Require_input_range std::make_pair

typedef double Real;
using namespace std;

class IInput_parser {
    public:
        typedef vector<string> tokenized_line;

        IInput_parser() = default;
        ~IInput_parser() = default;

        enum ERROR {
            UNSUPPORTED_DATATYPE,
            INCORRECT_TOKEN_COUNT,
        };

        virtual vector<tokenized_line> parse() = 0;
};

class Input_parser : public IInput_parser {

  public:

    Input_parser(string filename);
    vector<tokenized_line> parse();
    void split(const string& line, const char& delimiter, tokenized_line& target);

  private:

    vector<tokenized_line> m_settings;
    std::deque<string> m_raw_file;
    ifstream m_file;
    size_t m_num_starts{1};

    void read();
    bool is_not_start(const string& line);
    string discard_whitespace(const string& line);
    bool is_not_blank(const string& token);
    bool is_not_commented(const string& line);
    void assert_last_token_start(size_t line_no);
};

class Input {
  public:
  typedef vector<string> tokenized_line;
  typedef std::map<int, std::vector<tokenized_line>> Settings_per_start;

  typedef map<string,string> parameter_value_pair;
  typedef map<string, parameter_value_pair> brand_settings_map;
  typedef map<string, brand_settings_map> key_settings_map;
  typedef map<int, key_settings_map> start_settings_map;

  vector<string> KEYS;
	vector<string> SysList;
	vector<string> MolList;
	vector<string> MonList;
	vector<string> LatList;
	vector<string> AliasList;
	vector<string> NewtonList;
	vector<string> OutputList;
	vector<string> MesodynList;
	vector<string> ClengList;
	vector<string> TengList;
	vector<string> VarList;
	vector<string> StateList;
	vector<string> ReactionList;
  IInput_parser* m_parser;
  std::vector<tokenized_line> m_file_content;
  string name;
  map<string,string> parameter_value_pairs;
  start_settings_map m_settings;
  std::map<const std::string, vector<std::string>*> ListMap;
  std::map<const std::string, std::pair<int, int>> RangeMap;
  vector<string> out_options;

  Input(std::string filename);
  ~Input() { }

  void PreProcess();

  OutputInfo output_info;

  enum class INPUT_ERROR {
    ALREADY_DEFINED,
    UNKNOWN,
    OUT_NAME_NOT_FOUND,
    INPUT_RANGE,
    OUTPUT_FOLDER,
  };

  enum TOKENS {
      LINE,
      KEY,
      BRAND,
      PARAMETER,
      VALUE,
  };
  
  template<typename T>
  bool InSet(vector<T>& KEYS, T key) {
    auto iter = find(KEYS.begin(), KEYS.end(), key);
    return iter != KEYS.end();
  }

	template<typename T>
	bool InSet(vector<T> & KEYS, int &pos, T key) {
    auto iter = find(KEYS.begin(), KEYS.end(), key);
    pos = (int)std::distance(KEYS.begin(), iter);
    return iter != KEYS.end();
  }


  //Brackets verification for molecule compositions
  bool EvenDelimiters(string exp, vector<int>& open, vector<int>& close, char delimiter_open, char delimiter_close);
  bool EvenSquareBrackets(string exp, vector<int> &open, vector<int> &close);
  bool EvenBrackets(string exp, vector<int> &open, vector<int> &close);

  //Input validation
  bool ValidateKey(const std::vector<std::string>& KEYS, const std::string& key);
  void AssertRequestedOutputKeyExists(int start, const std::string& requested_key);
  bool OutputTypeRequested(int start, const string& output_option);
  void AssertRequestedOutputBrandExists(int start, const string& requested_key, const string& requested_brand);
  bool CheckInput(void);
  string LineError(const tokenized_line&);
  string PrintList(const std::vector<std::string>&);
  int GetNumStarts();
  void ValidateBrands();
  void ValidateInputKeys();
  bool CheckParameters(string key, string brand, int start, std::vector<std::string> &KEYS, std::vector<std::string> &PARAMETERS,std::vector<std::string> &VALUES);


  //List building
  bool MakeLists(int start);
  bool VerifyRange(int, int, int, std::string&);
  void ForceNonEmptyNameFor(const std::string, const std::string...);
  void ForceNonEmptyNameFor(std::string);
  bool LoadList(std::vector<std::string> &, const string,int, int, int);

  //Loading strings and converting to variables
  bool LoadItems(int start, string template_,std::vector<std::string> &Out_key, std::vector<std::string> &Out_name, std::vector<std::string> &Out_prop);
  void split(const string& line, const char& delimiter, tokenized_line& target);
  bool IsDigit( string token );

  template<typename Datatype>
  Datatype to_value(const string& PARAMETER, Datatype& target) {
        Datatype value;
        std::istringstream buffer{PARAMETER};
        buffer >> value;
        target = value;
        return value;
  }

  bool to_value(const string& PARAMETER, bool& target) {
        bool value;
        std::istringstream buffer{PARAMETER};
        buffer >> std::boolalpha >> value;
        target = value;
        return value;
  }

  template<typename Datatype>
  Datatype to_value(const string& PARAMETER) {
        Datatype value;
        std::istringstream buffer{PARAMETER};
        buffer >> value;
        return value;
  }

  bool to_value(const string& PARAMETER) {
        bool value;
        std::istringstream buffer{PARAMETER};
        buffer >> std::boolalpha >> value;
        return value;
  }

  bool ReadFile(string fname, string &In_buffer);

  // Legacy code used for compatability, please use the templates above.
  bool is_bool( string token );
  int Get_int(string, int );
	bool Get_int(string, int &, const std::string &);
	bool Get_int(string, int &, int, int, const std::string &);
	string Get_string(string, const string &);
	bool Get_string(string, string &, const  std::string & );
	bool Get_string(string , string &, std::vector<std::string>&, const std::string &);
	Real Get_Real(string, Real );
	bool Get_Real(string, Real &, const std::string &);
	bool Get_Real(string, Real &, Real, Real, const std::string &);
	bool Get_bool(string, bool );
	bool Get_bool(string, bool &, const std::string &);

  private:
	void parseOutputInfo();

};

#endif