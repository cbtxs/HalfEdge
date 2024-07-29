#include <ranges>
#include <vector>
#include <iostream>

template <typename Rng>
struct my_view : std::ranges::view_base {
    Rng base_;

    my_view() = default;
    my_view(Rng base) : base_(std::move(base)) {}

    auto begin() { return base_.begin(); }
    auto end() { return base_.end(); }
};

int main() {
    std::vector<int> vec = {1, 2, 3, 4, 5};
    my_view view(vec);

    for (int i : view) {
        std::cout << i << " ";
    }
}

