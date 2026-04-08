#include <pybind11/pybind11.h>
#include <string>

int run_contamination(
    int G, int K, int N_nei, int N, int N_MB, int N_tail, int n_record, int seed,
    const std::string& output_list,
    const std::string& data_name,
    const std::string& nei_name,
    const std::string& dist_name,
    const std::string& label_name,
    const std::string& cell_size_name,
    const std::string& MB_dir
);

namespace py = pybind11;

PYBIND11_MODULE(deContamination, m) {
    m.doc() = "deContamination Python bindings";

    m.def(
        "run",
        &run_contamination,
        py::arg("G"),
        py::arg("K"),
        py::arg("N_nei"),
        py::arg("N"),
        py::arg("N_MB"),
        py::arg("N_tail"),
        py::arg("n_record"),
        py::arg("seed"),
        py::arg("output_list"),
        py::arg("data_name"),
        py::arg("nei_name"),
        py::arg("dist_name"),
        py::arg("label_name"),
        py::arg("cell_size_name"),
        py::arg("MB_dir")
    );
}