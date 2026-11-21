// Matched outer covers on one cyclic, breadth-first material-box hierarchy.
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Integer = std::int64_t;
using Cell = std::uint64_t;
constexpr int kMaterial = 6;
constexpr int kDirect = 729; // 3^6
constexpr int kLinear = 64; // 2^6
constexpr int kReadingParts = 8; // six real residual and two denominator parts
constexpr int kCircuit = 81; // 3^4
constexpr int kCircuits = 15;
constexpr int kPower[] = {1, 3, 9, 27, 81, 243};

template <int N>
using Net = std::array<Integer, N>;

enum class Status { unresolved, outside, inside };

Integer floor_half(Integer x) {
    return x / 2 - (x < 0 && x % 2 != 0);
}

template <int N>
void split(const Net<N>& net, Net<N>* low, Net<N>* high, int axis) {
    const int stride = kPower[axis];
    for (int i = 0; i < N; i += 3 * stride) {
        for (int j = i; j < i + stride; ++j) {
            const auto a = net[j], b = net[j + stride], c = net[j + 2 * stride];
            const auto ab = floor_half(a + b), bc = floor_half(b + c);
            const auto midpoint = floor_half(ab + bc);
            if (low) {
                (*low)[j] = a;
                (*low)[j + stride] = ab;
                (*low)[j + 2 * stride] = midpoint;
            }
            if (high) {
                (*high)[j] = midpoint;
                (*high)[j + stride] = bc;
                (*high)[j + 2 * stride] = c;
            }
        }
    }
}

bool positive_quadratic(Integer a, Integer b, Integer c) {
    // Controls are at most 55 bits. Products fit in signed 128 bits.
    return a > 0 && c > 0 && (b >= 0 || (__int128)a * c > (__int128)b * b);
}

template <int N>
bool whole_coordinate_positive(const Net<N>& net, int axis) {
    const int stride = kPower[axis];
    for (int i = 0; i < N; i += 3 * stride) {
        for (int j = i; j < i + stride; ++j) {
            if (!positive_quadratic(net[j], net[j + stride], net[j + 2 * stride])) {
                return false;
            }
        }
    }
    return true;
}

using LinearNet = std::array<Net<kLinear>, kReadingParts>;

void split_linear(const LinearNet& net, LinearNet& low, LinearNet& high, int axis) {
    const int stride = 1 << (kMaterial - 1 - axis);
    for (int component = 0; component < kReadingParts; ++component) {
        for (int i = 0; i < kLinear; i += 2 * stride) {
            for (int j = i; j < i + stride; ++j) {
                const auto a = net[component][j], b = net[component][j + stride];
                const auto midpoint = floor_half(a + b);
                low[component][j] = a;
                low[component][j + stride] = midpoint;
                high[component][j] = midpoint;
                high[component][j + stride] = b;
            }
        }
    }
}

std::array<std::array<double, 2 * kMaterial>, kMaterial> centered_preconditioner(
    const LinearNet& lower, const LinearNet& negative_upper) {
    std::array<std::array<double, 2 * kMaterial>, kMaterial> matrix{};
    for (int component = 0; component < kMaterial; ++component) {
        for (int corner = 0; corner < kLinear; ++corner) {
            // A binary corner moment gives a center derivative up to one common factor.
            const auto value = lower[component][corner] - negative_upper[component][corner];
            for (int axis = 0; axis < kMaterial; ++axis) {
                const int sign = ((corner >> (kMaterial - 1 - axis)) & 1) ? 1 : -1;
                matrix[component][axis] += sign * static_cast<double>(value);
            }
        }
        matrix[component][component + kMaterial] = 1.;
    }
    for (int j = 0; j < kMaterial; ++j) {
        const auto pivot = std::max_element(
            matrix.begin() + j, matrix.end(),
            [j](const auto& a, const auto& b) { return std::abs(a[j]) < std::abs(b[j]); }
        );
        const double scale = (*pivot)[j];
        if (scale == 0. || !std::isfinite(scale)) return {};
        matrix[j].swap(*pivot);
        for (double& x : matrix[j]) x /= scale;
        for (int i = 0; i < kMaterial; ++i) {
            if (i == j) continue;
            const double factor = matrix[i][j];
            for (int k = 0; k < 2 * kMaterial; ++k) {
                matrix[i][k] -= factor * matrix[j][k];
            }
        }
    }
    return matrix;
}

__int128 square_integer(Integer x) {
    return (__int128)x * x;
}

bool centered_excludes(const LinearNet& lower, const LinearNet& upper) {
    __int128 residual_lower2 = 0, denominator_upper2 = 0;
    for (int component = 0; component < kReadingParts; ++component) {
        const auto [a, b] = std::minmax_element(lower[component].begin(), lower[component].end());
        const auto [c, d] = std::minmax_element(upper[component].begin(), upper[component].end());
        if (component < kMaterial) {
            residual_lower2 += square_integer(std::max({Integer(0), *a, *c}));
        } else {
            denominator_upper2 += square_integer(std::max(-*a, -*c));
        }
    }
    if (residual_lower2 > denominator_upper2) return true;

    const auto preconditioner = centered_preconditioner(lower, upper);
    for (const auto& row : preconditioner) {
        double norm = 0.;
        for (int j = 0; j < kMaterial; ++j) norm = std::max(norm, std::abs(row[kMaterial + j]));
        if (norm == 0. || !std::all_of(row.begin() + kMaterial, row.end(),
                                      [](double x) { return std::isfinite(x); })) continue;
        std::array<Integer, kMaterial> direction;
        Integer direction_norm2 = 0;
        for (int j = 0; j < kMaterial; ++j) {
            direction[j] = std::lround(row[kMaterial + j] / norm * (1 << 16));
            direction_norm2 += direction[j] * direction[j];
        }
        Integer plus = INT64_MAX, minus = INT64_MAX;
        for (int k = 0; k < kLinear; ++k) {
            Integer lo = 0, neg_hi = 0;
            for (int j = 0; j < kMaterial; ++j) {
                const auto q = direction[j];
                lo += (q >= 0) ? q * lower[j][k] : -q * upper[j][k];
                neg_hi += (q >= 0) ? q * upper[j][k] : -q * lower[j][k];
            }
            plus = std::min(plus, lo);
            minus = std::min(minus, neg_hi);
        }
        const auto distance = std::max({Integer(0), plus, minus});
        if (square_integer(distance) > direction_norm2 * denominator_upper2) return true;
    }
    return false;
}

struct DirectNode {
    Cell path;
    Net<kDirect> outside;
    Net<kDirect> inside;
};

struct Direct {
    using Node = DirectNode;
    Node root;
    std::uint64_t polynomial_tests = 0;

    explicit Direct(std::istream& in) {
        root.path = 1;
        in.read(reinterpret_cast<char*>(&root.outside), sizeof(root.outside));
        in.read(reinterpret_cast<char*>(&root.inside), sizeof(root.inside));
    }

    Status test(const Node& node, int /* depth */, int next_axis) {
        polynomial_tests += 2;
        // Let the direct baseline use the same quadratic positivity lemma along
        // the next axis, in addition to its complete joint noise score.
        const int axis = kMaterial - 1 - next_axis;
        if (whole_coordinate_positive<kDirect>(node.outside, axis)) return Status::outside;
        if (whole_coordinate_positive<kDirect>(node.inside, axis)) return Status::inside;
        return Status::unresolved;
    }

    void push(const Node& node, std::vector<Node>& target, int split_axis) {
        const auto i = target.size();
        target.resize(i + 2);
        auto& low = target[i];
        auto& high = target[i + 1];
        low.path = node.path * 2;
        high.path = low.path + 1;
        const int axis = kMaterial - 1 - split_axis; // NumPy C order.
        split<kDirect>(node.outside, &low.outside, &high.outside, axis);
        split<kDirect>(node.inside, &low.inside, &high.inside, axis);
    }

    void end_level(int) {}
};

struct ChartNode { Cell path; };

struct CenteredNode {
    DirectNode score;
    LinearNet lower;
    LinearNet negative_upper;
    Cell path;
};

struct Centered {
    using Node = CenteredNode;
    Direct direct;
    Node root;
    std::uint64_t polynomial_tests = 0;
    std::uint64_t centered_tests = 0;
    std::uint64_t centered_exclusions = 0;

    explicit Centered(std::istream& in) : direct(in) {
        root.path = 1;
        root.score = direct.root;
        in.read(reinterpret_cast<char*>(&root.lower), sizeof(root.lower));
        in.read(reinterpret_cast<char*>(&root.negative_upper), sizeof(root.negative_upper));
    }

    Status test(const Node& node, int depth, int next_axis) {
        const auto score = direct.test(node.score, depth, next_axis);
        polynomial_tests = direct.polynomial_tests;
        if (score != Status::unresolved) return score;
        ++centered_tests;
        if (centered_excludes(node.lower, node.negative_upper)) {
            ++centered_exclusions;
            return Status::outside;
        }
        return Status::unresolved;
    }

    void push(const Node& node, std::vector<Node>& target, int axis) {
        const auto i = target.size();
        target.resize(i + 2);
        auto& low = target[i];
        auto& high = target[i + 1];
        low.path = node.path * 2;
        high.path = low.path + 1;
        low.score.path = low.path;
        high.score.path = high.path;
        split<kDirect>(node.score.outside, &low.score.outside, &high.score.outside, kMaterial - 1 - axis);
        split<kDirect>(node.score.inside, &low.score.inside, &high.score.inside, kMaterial - 1 - axis);
        split_linear(node.lower, low.lower, high.lower, axis);
        split_linear(node.negative_upper, low.negative_upper, high.negative_upper, axis);
    }

    void end_level(int) {}
};

struct Projection {
    std::array<int, kMaterial> stride_axis;
    Net<kCircuit> root;

    struct Record {
        Net<kCircuit> net;
        bool excluded;
    };

    std::unordered_map<Cell, Record> current;
    std::unordered_map<Cell, Record> previous;

    Projection() = default;

    void read(std::istream& in) {
        std::array<std::int8_t, kMaterial> axes;
        in.read(reinterpret_cast<char*>(&axes), sizeof(axes));
        std::copy(axes.begin(), axes.end(), stride_axis.begin());
        in.read(reinterpret_cast<char*>(&root), sizeof(root));
    }

    void advance(int axis) {
        if (stride_axis[axis] >= 0) {
            previous = std::move(current);
            current.clear();
        }
    }

    Cell projected_path(Cell path, int depth, const std::array<int, kMaterial>& order) const {
        Cell out = 1;
        for (int level = 0; level < depth; ++level) {
            if (stride_axis[order[level % kMaterial]] >= 0) {
                out = (out << 1) | ((path >> (depth - 1 - level)) & 1);
            }
        }
        return out;
    }

    bool excludes(Cell path, int depth, const std::array<int, kMaterial>& order,
                  std::uint64_t& polynomial_tests) {
        const auto key = projected_path(path, depth, order);
        const auto exists = current.find(key);
        if (exists != current.end()) return exists->second.excluded;

        Record value{};
        const auto parent = previous.find(key / 2);
        if (parent != previous.end()) {
            const auto axis = [&] {
                int level = depth - 1;
                while (stride_axis[order[level % kMaterial]] < 0) --level;
                return stride_axis[order[level % kMaterial]];
            }();
            if (key % 2 == 0) split<kCircuit>(parent->second.net, &value.net, nullptr, axis);
            else split<kCircuit>(parent->second.net, nullptr, &value.net, axis);
        } else {
            value.net = root;
            for (int level = 0; level < depth; ++level) {
                const int axis = stride_axis[order[level % kMaterial]];
                if (axis < 0) continue;
                Net<kCircuit> next;
                if (((path >> (depth - 1 - level)) & 1) == 0) {
                    split<kCircuit>(value.net, &next, nullptr, axis);
                } else {
                    split<kCircuit>(value.net, nullptr, &next, axis);
                }
                value.net = next;
            }
        }
        // A positive degree-two slice in any of its four axes excludes
        // every value of that coordinate and both other physical coordinates.
        for (int axis = 0; axis < 4 && !value.excluded; ++axis) {
            ++polynomial_tests;
            value.excluded = whole_coordinate_positive<kCircuit>(value.net, axis);
        }
        current.emplace(key, value);
        return value.excluded;
    }
};

struct Chart {
    using Node = ChartNode;
    ChartNode root{1};
    std::array<Projection, kCircuits> projection;
    const std::array<int, kMaterial>& order;
    std::uint64_t polynomial_tests = 0;

    Chart(std::istream& in, const std::array<int, kMaterial>& axis_order) : order(axis_order) {
        for (auto& p : projection) p.read(in);
    }

    Status test(const Node& node, int depth, int /* next axis */) {
        for (auto& p : projection) {
            if (p.excludes(node.path, depth, order, polynomial_tests)) return Status::outside;
        }
        return Status::unresolved;
    }

    void push(const Node& node, std::vector<Node>& target, int /* split axis */) {
        target.push_back({2 * node.path});
        target.push_back({2 * node.path + 1});
    }

    void end_level(int axis) {
        for (auto& p : projection) p.advance(axis);
    }
};

struct Checkpoint {
    std::uint64_t visited = 0;
    std::uint64_t omitted = 0;
    std::uint64_t inside = 0;
    std::uint64_t retained = 1;
    double outer_volume = 1.;
    double inner_volume = 0.;
};

template <class Method>
void cover(Method& method, const std::array<int, kMaterial>& axis_order,
           std::uint64_t max_visits, int max_depth, const std::string& output) {
    const auto wall_start = std::chrono::steady_clock::now();
    const auto cpu_start = std::clock();
    Checkpoint state;
    std::uint64_t previous = ~std::uint64_t(0);
    const std::array<double, 29> thresholds{
        .99, .9, .5, .3, .1, .03, .01, .003, .001, .0003, .0001,
        .00003, .00001, .000003, .000001, .0000003, .0000001,
        3.e-8, 1.e-8, 3.e-9, 1.e-9, 3.e-10, 1.e-10, 3.e-11, 1.e-11,
        3.e-12, 1.e-12, 3.e-13, 1.e-13
    };
    size_t threshold = 0;
    std::vector<typename Method::Node> current{method.root}, next;
    std::vector<Cell> inside_cells, terminal;
    int depth = 0;
    std::size_t i = 0;

    auto snapshot = [&] {
        if (previous == state.visited) return;
        previous = state.visited;
        const double cpu = double(std::clock() - cpu_start) / CLOCKS_PER_SEC;
        const double wall =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
        std::cout << std::setprecision(17)
                  << "{\"visited\":" << state.visited
                  << ",\"excluded\":" << state.omitted
                  << ",\"inside\":" << state.inside
                  << ",\"retained\":" << state.retained
                  << ",\"undecided\":" << state.retained - state.inside
                  << ",\"outer_fraction\":" << state.outer_volume
                  << ",\"inner_fraction\":" << state.inner_volume
                  << ",\"depth\":" << depth
                  << ",\"polynomial_tests\":" << method.polynomial_tests
                  << ",\"cpu_seconds_sieve\":" << cpu
                  << ",\"wall_seconds_sieve\":" << wall;
        if constexpr (requires { method.centered_tests; }) {
            std::cout << ",\"centered_tests\":" << method.centered_tests
                      << ",\"centered_exclusions\":" << method.centered_exclusions;
        }
        std::cout << "}\n";
    };

    snapshot();
    for (; !current.empty() && state.visited < max_visits; ++depth) {
        const int axis = axis_order[depth % kMaterial];
        const double volume = std::ldexp(1., -depth);
        i = 0;
        for (; i < current.size() && state.visited < max_visits; ++i) {
            const auto& node = current[i];
            ++state.visited;
            const auto verdict = method.test(node, depth, axis);
            if (verdict == Status::outside) {
                --state.retained;
                ++state.omitted;
                state.outer_volume -= volume;
            } else if (verdict == Status::inside) {
                ++state.inside;
                state.inner_volume += volume;
                inside_cells.push_back(node.path);
            } else if (depth < max_depth) {
                ++state.retained;
                method.push(node, next, axis);
            } else {
                terminal.push_back(node.path);
            }
            bool quality_event = false;
            while (threshold < thresholds.size() && state.outer_volume <= thresholds[threshold]) {
                ++threshold;
                quality_event = true;
            }
            if (quality_event || state.visited == 1 || state.visited == 100 ||
                state.visited % 1000 == 0 || state.visited == max_visits) {
                snapshot();
            }
        }
        snapshot();
        if (i < current.size()) break;
        method.end_level(axis);
        current.swap(next);
        next.clear();
        i = 0;
    }

    const auto count = terminal.size() + current.size() - i + next.size();
    if (count + inside_cells.size() != state.retained) {
        throw std::logic_error("Frontier leaf accounting failed");
    }
    std::ofstream file(output, std::ios::binary);
    const std::array<std::uint64_t, 2> lengths{inside_cells.size(), count};
    file.write(reinterpret_cast<const char*>(&lengths), sizeof(lengths));
    file.write(reinterpret_cast<const char*>(inside_cells.data()), inside_cells.size() * sizeof(Cell));
    file.write(reinterpret_cast<const char*>(terminal.data()), terminal.size() * sizeof(Cell));
    for (; i < current.size(); ++i) {
        file.write(reinterpret_cast<const char*>(&current[i].path), sizeof(Cell));
    }
    for (const auto& node : next) file.write(reinterpret_cast<const char*>(&node.path), sizeof(Cell));
    if (!file) throw std::runtime_error("Could not write cell archive");
}

int main(int argc, char** argv) {
    if (argc != 6) throw std::invalid_argument("arm input output max_visits max_depth");
    std::ifstream in(argv[2], std::ios::binary);
    std::array<int, kMaterial> order;
    for (auto& axis : order) axis = in.get();
    auto check = order;
    std::sort(check.begin(), check.end());
    if (check != std::array{0, 1, 2, 3, 4, 5}) throw std::runtime_error("Invalid axis order");
    const auto max_visits = std::stoull(argv[4]);
    const int max_depth = std::stoi(argv[5]);
    if (max_depth < 0 || max_depth > 48) throw std::out_of_range("Depth outside [0,48]");
    const std::string arm = argv[1];
    if (arm == "direct") {
        Direct method(in);
        if (!in) throw std::invalid_argument("Missing direct coefficients");
        cover(method, order, max_visits, max_depth, argv[3]);
    } else if (arm == "direct_centered") {
        Centered method(in);
        if (!in) throw std::invalid_argument("Missing centered coefficients");
        cover(method, order, max_visits, max_depth, argv[3]);
    } else if (arm == "chart") {
        Chart method(in, order);
        if (!in) throw std::invalid_argument("Missing chart coefficients");
        cover(method, order, max_visits, max_depth, argv[3]);
    } else {
        throw std::invalid_argument("Unknown cover arm");
    }
}
