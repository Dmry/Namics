#include "file_writer.h"

Register_class<IProfile_writer, Vtk_structured_grid_writer, Writable_filetype, Lattice*, Writable_file> Vtk_structured_grid_writer_factory(Writable_filetype::VTK_STRUCTURED_GRID);
Register_class<IProfile_writer, Vtk_structured_points_writer, Writable_filetype, Lattice*, Writable_file> Vtk_structured_points_writer_factory(Writable_filetype::VTK_STRUCTURED_POINTS);
Register_class<IProfile_writer, Pro_writer, Writable_filetype, Lattice*, Writable_file> Pro_writer_factory(Writable_filetype::PRO);

map<std::string, Writable_filetype> Profile_writer::input_options {
        {"vtk_structured_grid", Writable_filetype::VTK_STRUCTURED_GRID},
        {"vtk_structured_points", Writable_filetype::VTK_STRUCTURED_POINTS},
        {"pro", Writable_filetype::PRO}
    };

std::map<Writable_filetype, std::string> Writable_file::extension_map {
    {Writable_filetype::KAL, "kal"},
    {Writable_filetype::VTK_STRUCTURED_GRID, "vtk"},
    {Writable_filetype::VTK_STRUCTURED_POINTS, "vtk"},
    {Writable_filetype::CSV, "csv"},
    {Writable_filetype::PRO, "pro"},
    {Writable_filetype::JSON, "json"}
};

Writable_file::Writable_file(const std::string filename_, Writable_filetype filetype_, int identifier_)
: m_identifier{identifier_}, m_filename{filename_}, m_filetype{filetype_}
{
    append_extension();
}

Writable_file::~Writable_file()
{

}

void Writable_file::increment_identifier()
{
    string dot = "\\.";
    string s_regex = "(\\_" + to_string(m_identifier) + dot + Writable_file::extension_map[m_filetype] + ')';
    std::regex ex(s_regex);

    string replacement = "_"+ to_string(++m_identifier) + "." + Writable_file::extension_map[m_filetype];

    m_filename = regex_replace(m_filename,ex, replacement);
}

Writable_filetype Writable_file::get_filetype()
{
    return m_filetype;
}

string Writable_file::get_filename()
{
    return m_filename;
}

void Writable_file::append_extension()
{
    m_filename.append("_" + to_string(m_identifier) + "." + Writable_file::extension_map[m_filetype]);
}

IProfile_writer::IProfile_writer(Lattice* geometry_, Writable_file file_)
: m_geometry{geometry_}, m_file{file_}, m_adapter{geometry_}
{
    m_filestream.precision(configuration.precision);
}

IProfile_writer::~IProfile_writer()
{

}

void IProfile_writer::bind_data(map< string, shared_ptr<IOutput_ptr>>& profiles_)
{
    m_profiles.insert(profiles_.begin(), profiles_.end());
}

void IProfile_writer::bind_subystem_loop(Boundary_mode mode_) {
    switch (mode_) {
        case Boundary_mode::WITH_BOUNDS:
            subsystem_loop = std::bind(&Lattice_accessor::system_plus_bounds, m_adapter, std::placeholders::_1);
            break;
        case Boundary_mode::WITHOUT_BOUNDS:
            subsystem_loop = std::bind(&Lattice_accessor::skip_bounds, m_adapter, std::placeholders::_1);
            break;
    }
}


Vtk_structured_grid_writer::Vtk_structured_grid_writer(Lattice* geometry_, Writable_file file_)
: IProfile_writer(geometry_, file_)
{
    configuration.boundary_mode = Boundary_mode::WITHOUT_BOUNDS;
}

Vtk_structured_grid_writer::~Vtk_structured_grid_writer()
{

}

void Vtk_structured_grid_writer::prepare_for_data()
{
    bind_subystem_loop(configuration.boundary_mode);

    m_filestream.open(m_file.get_filename(), std::ios_base::out);

	std::ostringstream vtk;

    int MX = m_geometry->MX;
	int MY = m_geometry->MY;
	int MZ = m_geometry->MZ;

	vtk << "# vtk DataFile Version 4.2 \n";
	vtk << "VTK output \n";
	vtk << "ASCII\n";
	vtk << "DATASET STRUCTURED_GRID \n";
	vtk << "DIMENSIONS " << MX << " " << MY << " " << MZ << "\n";
	vtk << "POINTS " << MX * MY * MZ << " int\n";

	for (int x = 1; x < MX + 1; ++x)
		for (int y = 1; y < MY + 1; ++y)
			for (int z = 1 ; z < MZ + 1 ; ++z )
				vtk << x << " " << y << " " << z << "\n";

	vtk << "POINT_DATA " << MX * MY * MZ << "\n";

	m_filestream << vtk.str();
	m_filestream.flush();

	m_filestream.close();
}

void Vtk_structured_grid_writer::write()
{
	m_filestream.open(m_file.get_filename(), std::ios_base::app);

    for (auto& profile : m_profiles) {

	    std::ostringstream vtk;

	    vtk << "SCALARS " << profile.first << " float\nLOOKUP_TABLE default \n";

	    subsystem_loop(
            [this, &vtk, profile] (size_t x, size_t y, size_t z)
            {
	    		vtk << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
            }
        );

	    m_filestream << vtk.str();
	    m_filestream.flush();

    }

	m_filestream.close();
    m_file.increment_identifier();
}

Vtk_structured_points_writer::Vtk_structured_points_writer(Lattice* geometry_, Writable_file file_)
: IProfile_writer(geometry_, file_)
{
    configuration.boundary_mode = Boundary_mode::WITHOUT_BOUNDS;
}

Vtk_structured_points_writer::~Vtk_structured_points_writer()
{

}

void Vtk_structured_points_writer::prepare_for_data()
{
    bind_subystem_loop(configuration.boundary_mode);

    m_filestream.open(m_file.get_filename(), std::ios_base::out);

	std::ostringstream vtk;

    int MX = m_geometry->MX;
	int MY = m_geometry->MY;
	int MZ = m_geometry->MZ;

    vtk << "# vtk DataFile Version 4.2 \n";
	vtk << "VTK output \n";
	vtk << "ASCII\n";
	vtk << "DATASET STRUCTURED_POINTS \n";
	vtk << "DIMENSIONS " << MX << " " << MY << " " << MZ << "\n";

    m_filestream << vtk.str();
	m_filestream.flush();

	m_filestream.close();
}

void Vtk_structured_points_writer::write()
{
    m_filestream.open(m_file.get_filename(), std::ios_base::app);

    for (auto& profile : m_profiles) {

	    std::ostringstream vtk;

        vtk << "SCALARS " << profile.first << " float\nLOOKUP_TABLE default \n";

        #ifdef PAR_MESODYN
            stl::host_vector<Real> m_data;

            subsystem_loop(
                [this, &vtk, profile] (size_t x, size_t y, size_t z) {
	    		    vtk << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
                }
            );
        #else
        subsystem_loop(
            [this, &vtk, profile] (size_t x, size_t y, size_t z) {
	    		vtk << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
            }
        );
        #endif

        
        m_filestream << vtk.str();
	    m_filestream.flush();
    }

	m_filestream.close();
    m_file.increment_identifier();
}

Pro_writer::Pro_writer(Lattice* geometry_, Writable_file file_)
: IProfile_writer(geometry_, file_)
{
    configuration.boundary_mode = Boundary_mode::WITH_BOUNDS;
}

Pro_writer::~Pro_writer()
{

}

void Pro_writer::prepare_for_data()
{
    bind_subystem_loop(configuration.boundary_mode);
    
    m_filestream.open(m_file.get_filename(), std::ios_base::out);

	std::ostringstream pro;

    pro << "x";

    if (m_geometry->gradients > 1)
    {
        pro << "\ty";
    }   

    if (m_geometry->gradients > 2)
    {
        pro << "\tz";
    }

    pro << '\n';

    subsystem_loop(
        [this, &pro] (size_t x, size_t y, size_t z) {
			pro << x;
            if (m_geometry->gradients > 1)
            {
                pro << "\t" << y;
                if (m_geometry->gradients > 1)
                {
                    pro << "\t" << z;
                }

            }
            pro << '\n';
        }
    );

    m_filestream << pro.str();
	m_filestream.flush();

	m_filestream.close();
}

void Pro_writer::write()
{
    string temp_filename = m_file.get_filename()+".temp";

    for (auto& profile : m_profiles) {
        m_filestream.open(temp_filename, std::ios_base::out);
        std::ifstream in_file(m_file.get_filename());
        if (!in_file) {
            std::cerr << "Could not open input file. Did you prepare for data?\n";
            throw 1;
        }

	    std::ostringstream pro;

        std::string in;
        std::getline(in_file, in);
        m_filestream << in << "\t" << profile.first << '\n';

        subsystem_loop(
            [this, &pro, &in, &in_file, profile] (size_t x, size_t y, size_t z) {
                std::getline(in_file, in);
                m_filestream << in << "\t" << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
                m_filestream << pro.str();
	            m_filestream.flush();
            }
        );

        in_file.close();
        m_filestream.close();

        rename(temp_filename.c_str(), m_file.get_filename().c_str());
    }

	m_filestream.close();
    m_file.increment_identifier();
}

IParameter_writer::IParameter_writer(Writable_file file_)
: m_file{file_}
{
    m_filestream.precision(configuration.precision);
}

void IParameter_writer::bind_data(map< string, shared_ptr<IOutput_ptr>>& params_)
{
    m_params.insert(params_.begin(), params_.end());
}

void IParameter_writer::register_categories(std::map<string, CATEGORY>& category_) {
    m_categories = category_;
}

Kal_writer::Kal_writer(Writable_file file_)
: IParameter_writer(file_)
{
    
}

Kal_writer::~Kal_writer() {

}

void Kal_writer::write() {
    m_filestream.open(m_file.get_filename(), ios_base::app);
    for (auto& description : m_selected_variables) {
        m_filestream << m_params[description]->data() << '\t';
    }
    m_filestream << '\n';
    m_filestream.close();
}

void Kal_writer::prepare_for_data(vector<string>& selected_variables_) {
    m_selected_variables = selected_variables_;

    m_filestream.open(m_file.get_filename());

    for (auto& description : m_selected_variables)
        m_filestream << description << '\t';
    m_filestream << '\n';

    m_filestream.close();
}

Csv_parameter_writer::Csv_parameter_writer(Writable_file file_)
: IParameter_writer(file_)
{
    
}

Csv_parameter_writer::~Csv_parameter_writer() {

}

void Csv_parameter_writer::write() {
    m_filestream.open(m_file.get_filename(), ios_base::app);
    for (auto& description : m_selected_variables) {
        m_filestream << m_params[description]->data() << ',';
    }
    m_filestream << '\n';
    m_filestream.close();
}

void Csv_parameter_writer::prepare_for_data(vector<string>& selected_variables_) {
    m_selected_variables = selected_variables_;

    m_filestream.open(m_file.get_filename());

    for (auto& description : m_selected_variables)
        m_filestream << description << ',';
    m_filestream << '\n';

    m_filestream.close();
}


JSON_parameter_writer::JSON_parameter_writer(Writable_file file_)
: IParameter_writer(file_), m_selected_constants(0), m_selected_metadata(0), m_selected_timespan(0)
{
    
}

JSON_parameter_writer::~JSON_parameter_writer() {
    if (m_state == STATE::IS_OPEN)
        finalize();
}

bool JSON_parameter_writer::is_number(std::string& value) {
    return std::regex_match(value, std::regex(R"(^[\d]{1,}(.[\d]{1,})?e?\-?[\d]+$)"));
}

void JSON_parameter_writer::prepare_for_data(std::vector<std::string>& selected_variables_) {
    m_selected_variables = selected_variables_;

    preprocess_categories();

    // Select variables from metadata and constants

    partition_selected_variables(m_metadata, m_selected_metadata);
    partition_selected_variables(m_constants, m_selected_constants);
    partition_selected_variables(m_timespan, m_selected_timespan);

    m_filestream.open(m_file.get_filename());

    open_json();

    write_list(m_selected_metadata);

    m_filestream << ',';

    open_object("Constants");
    write_list(m_selected_constants);
    close_object();

    open_array("Time Series");

    m_filestream << "{\n";
    write_list(m_selected_timespan);
    m_filestream << "}";

    m_filestream.close();
}

void JSON_parameter_writer::finalize() {
    m_filestream.open(m_file.get_filename(), ios_base::app);
    close_array();
    close_json();
    m_state = STATE::IS_CLOSED;
    m_filestream.close();
}

void JSON_parameter_writer::write() {
    m_filestream.open(m_file.get_filename(), ios_base::app);

    write_array_object(m_selected_timespan);

    m_filestream.close();
}

void JSON_parameter_writer::partition_selected_variables(std::vector<std::string>& partition, std::vector<std::string>& target) {
    auto it = std::partition(m_selected_variables.begin(), m_selected_variables.end(), [partition](std::string variable) {
        return std::find(partition.begin(), partition.end(), variable) != partition.end();
    });

    for ( auto itt = m_selected_variables.begin() ; itt != it ; ++itt)
        target.push_back(*itt);
}

void JSON_parameter_writer::preprocess_categories() {
    for (auto const& entry : m_categories)
        switch (entry.second) { //CATEGORY enum
            case IParameter_writer::CATEGORY::METADATA:
                m_metadata.push_back(entry.first);
                break;
            case IParameter_writer::CATEGORY::CONSTANT:
                m_constants.push_back(entry.first);
                break;
            case IParameter_writer::CATEGORY::TIMESPAN:
                m_timespan.push_back(entry.first);
                break;
        }
}

void JSON_parameter_writer::write_list(std::vector<std::string>& descriptions) {
    for (auto itt = descriptions.begin(); itt != descriptions.end(); ++itt) {
        m_filestream << "\"" << *itt << "\": ";

        string data =  m_params[*itt]->data();

        if (!is_number(data))
            m_filestream << "\"";
            
        m_filestream << m_params[*itt]->data();
        
        if (!is_number(data))
            m_filestream << "\"";

        if (itt != descriptions.end()-1)
            m_filestream << ",";

        m_filestream << "\n";
    }
}

JSON_parameter_writer::STATE JSON_parameter_writer::get_state() {
    return m_state;
}

void JSON_parameter_writer::write_array_object(std::vector<std::string>& descriptions) {
    m_filestream << ",\n{\n";
    write_list(descriptions);
    m_filestream << "}";
}

void JSON_parameter_writer::open_json() {
    m_state = STATE::IS_OPEN;
    m_filestream << "{\n";
}

void JSON_parameter_writer::close_json() {
    m_state = STATE::IS_CLOSED;
    m_filestream << "}\n";
}

void JSON_parameter_writer::open_object(const std::string& description) {
    m_filestream << "\"" << description << "\": {\n";
}

void JSON_parameter_writer::close_object() {
    m_filestream << "},\n";
}

void JSON_parameter_writer::open_array(const std::string& description) {
    m_filestream << "\"" << description << "\": [\n";
}

void JSON_parameter_writer::close_array() {
    m_filestream << "]\n";
}