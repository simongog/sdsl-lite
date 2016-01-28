#include "utils.hpp"
#include "index_types.hpp"
#include "collection.hpp"

#include "logging.hpp"

typedef struct cmdargs {
    std::string collection_dir;
} cmdargs_t;

void print_usage(const char* program)
{
    fprintf(stdout, "%s -c <collection dir>\n", program);
    fprintf(stdout, "where\n");
    fprintf(stdout, "  -c <collection dir>  : the collection dir.\n");
};

cmdargs_t parse_args(int argc, const char* argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    while ((op = getopt(argc, (char* const*)argv, "c:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
        }
    }
    if (args.collection_dir == "") {
        LOG(FATAL) << "Missing command line parameters.";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

template <class t_idx> void create_and_store(collection& col)
{
    using clock = std::chrono::high_resolution_clock;
    auto start = clock::now();
    t_idx idx(col);

    auto stop = clock::now();
    LOG(INFO) << "index construction in (s): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              / 1000.0f;
    auto output_file = col.path + "/index/index-" + idx.name() + "-" + sdsl::util::class_to_hash(idx) + ".sdsl";
    std::ofstream ofs(output_file);
    if (ofs.is_open()) {
        LOG(INFO) << "writing index to file : " << output_file;
        auto bytes = sdsl::serialize(idx, ofs);
        LOG(INFO) << "index size : " << bytes / (1024 * 1024) << " MiB";
        LOG(INFO) << "writing space usage visualization to file : " << output_file + ".html";
        std::ofstream vofs(output_file + ".html");
        sdsl::write_structure<sdsl::HTML_FORMAT>(vofs, idx);
    } else {
        LOG(FATAL) << "cannot write index to file : " << output_file;
    }
}

int main(int argc, const char* argv[])
{
    sdsl::construct_config::byte_algo_sa = sdsl::SE_SAIS;

    log::start_log(argc, argv, false);

    /* parse command line */
    cmdargs_t args = parse_args(argc, argv);

    /* parse collection directory */
    collection col(args.collection_dir);

    /* create index */
    {
        using index_type = INDEX_TYPE;
        create_and_store<index_type>(col);
    }

    return 0;
}
