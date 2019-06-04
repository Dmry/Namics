#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <algorithm>
#include <vector>
#include <regex>
#include <iostream>
#include <exception>
#include <map>
#include <functional>
#include <fstream>
#include "output_info.h"

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
    void split(const string& line, const string& delimiter, tokenized_line& target);

  private:

    vector<tokenized_line> m_settings;
    ifstream m_file;
    size_t m_num_starts{1};

    bool is_not_start(const string& line);
    string discard_whitespace(const string& line);
    bool is_not_blank(const string& token);
    bool is_not_commented(const string& line);
    bool assert_last_token_start(size_t line_no);
};

class Input {
  public:
    typedef vector<string> tokenized_line;

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

  Input(std::string filename);
  ~Input() { }

  OutputInfo output_info;

  enum TOKENS {
      LINE,
      KEY,
      BRAND,
      PARAMETER,
      VALUE,
  };

  IInput_parser* m_parser;
  const vector<tokenized_line> m_file_content;
  map<string,string> parameter_value_pair;
  
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

  bool CheckParameters(string key, string brand, int start, std::vector<std::string> &KEYS, std::vector<std::string> &PARAMETERS,std::vector<std::string> &VALUES);
  std::sregex_iterator NestedGroup(string& input, string& delimiter_open, string& delimiter_close);
  bool EvenSign(string input, string delimiter_open, string delimiter_close, vector<int>&, vector<int>&);
  bool EvenSquareBrackets(string exp, vector<int> &open, vector<int> &close);
  bool EvenBrackets(string exp, vector<int> &open, vector<int> &close);
  bool ValidateKeys();
  bool CheckInput(void);
  string name;
  bool MakeLists(int start);
  bool TestNum(std::vector<std::string> &S, string c,int num_low, int num_high, int UptoStartNumber );
  int GetNumStarts();
  bool IsDigit( string token );
  bool is_bool( string token );
  int Get_int(string, int );
	bool Get_int(string, int &, const std::string &);
  bool ArePair(char opening,char closing);
	bool Get_int(string, int &, int, int, const std::string &);
	string Get_string(string, const string &);
	bool Get_string(string, string &, const  std::string & );
	bool Get_string(string , string &, std::vector<std::string>&, const std::string &);
	Real Get_Real(string, Real );
	bool Get_Real(string, Real &, const std::string &);
	bool Get_Real(string, Real &, Real, Real, const std::string &);
	bool Get_bool(string, bool );
	bool Get_bool(string, bool &, const std::string &);
    void split(const string& line, const char& delimiter, tokenized_line& target);

  private:
	void parseOutputInfo();

};

#endif