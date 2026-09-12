/* Pure-C fixed cases and generated properties over borrowed kernel views. */
#include "duckvep_property.h"

GREATEST_MAIN_DEFS();

int main(int argc, char **argv) {
    GREATEST_MAIN_BEGIN();
#define DUCKVEP_PROPERTY_TEST(name) RUN_TEST(name);
#include "duckvep_property_tests.def"
#undef DUCKVEP_PROPERTY_TEST
    GREATEST_MAIN_END();
}
