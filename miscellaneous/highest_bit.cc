
int highest_bit(unsigned x) {
    return x == 0 ? -1 : 31 - __builtin_clz(x);
}


int highest_bit(uint64_t x) {
    return x == 0 ? -1 : 63 - __builtin_clzll(x);
}
