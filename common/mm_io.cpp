#include "mm_io.hpp"

Graph Read_MM_UD(const std::string& file_name)
{
    std::ifstream infile(file_name);
    if (!infile)
    {
        std::cerr << "Error: Cannot open file " << file_name << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    // コメント行をスキップし、サイズを読み取る
    do
    {
        getline(infile, line);
    } while (!infile.eof() && line[0] == '%');

    std::istringstream size_stream(line);
    int row_size, col_size, num_elements;
    size_stream >> row_size >> col_size >> num_elements;

    if (row_size != col_size)
    {
        std::cerr << "Error: The matrix is not square." << std::endl;
        exit(EXIT_FAILURE);
    }

    Graph G(row_size);
    WeightMap weightmap = get(boost::edge_weight, G);

    int row, col;
    double val;
    while (getline(infile, line))
    {
        if (line.empty() || line[0] == '%') continue;

        std::istringstream edge_stream(line);
        edge_stream >> row >> col >> val;

        // Matrix Market is 1-based index → convert to 0-based
        if (row > col)
        {
        auto e = add_edge(row - 1, col - 1, G).first;
        weightmap[e] = val;
        }
    }

    return G;
}

// ---- ユーティリティ：入力パスから拡張子を外したstemを得る
std::string file_stem(const std::string& path) {
    // 末尾の '/' or '\' の次からファイル名部分
    size_t pos = path.find_last_of("/\\");
    std::string name = (pos == std::string::npos) ? path : path.substr(pos+1);
    // 拡張子を落とす
    size_t dot = name.find_last_of('.');
    return (dot == std::string::npos) ? name : name.substr(0, dot);
}