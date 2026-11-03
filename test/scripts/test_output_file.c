/* Run from an empty directory inside the clone:
 * cc -std=c11 -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Werror \
 *    -I../../src/include ../scripts/test_output_file.c -o test_output_file
 * ./test_output_file
 */
#include "duckhts_output_file.h"

static int expect_contents(const char *path, const char *expected) {
    char buffer[128] = {0};
    FILE *file = fopen(path, "rb");
    if (!file) return -1;
    size_t length = fread(buffer, 1, sizeof(buffer) - 1, file);
    int failed = ferror(file);
    if (fclose(file) != 0) failed = 1;
    if (failed || length != strlen(expected) || strcmp(buffer, expected) != 0) return -1;
    return 0;
}

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "line %d: %s (errno %d)\n", __LINE__, #condition, errno); \
        return 1; \
    } \
} while (0)

int main(void) {
    duckhts_output_owner_t owner = {0};
    struct stat st;
    int fd;
    char *staged;
    umask(0077);

    fd = duckhts_output_open("output", 0, &owner);
    CHECK(fd >= 0);
    CHECK(write(fd, "original", 8) == 8);
    CHECK(close(fd) == 0);
    CHECK(duckhts_output_finish("output", &owner, 1) == 0);
    CHECK(stat("output", &st) == 0 && (st.st_mode & 0777) == 0600);
    CHECK(duckhts_output_open("output", 0, &owner) < 0 && errno == EEXIST);

    fd = duckhts_output_open("output", 1, &owner);
    CHECK(fd >= 0);
    staged = strdup(duckhts_output_path("output", &owner));
    CHECK(staged != NULL);
    CHECK(write(fd, "partial", 7) == 7);
    CHECK(close(fd) == 0);
    CHECK(expect_contents("output", "original") == 0);
    CHECK(duckhts_output_finish("output", &owner, 0) == 0);
    CHECK(expect_contents("output", "original") == 0);
    CHECK(access(staged, F_OK) < 0 && errno == ENOENT);
    free(staged);

    fd = duckhts_output_open("output", 1, &owner);
    CHECK(fd >= 0);
    CHECK(write(fd, "complete", 8) == 8);
    CHECK(close(fd) == 0);
    CHECK(duckhts_output_finish("output", &owner, 1) == 0);
    CHECK(expect_contents("output", "complete") == 0);
    CHECK(stat("output", &st) == 0 && (st.st_mode & 0777) == 0600);

    CHECK(mkdir("directory", 0700) == 0);
    fd = duckhts_output_open("directory", 1, &owner);
    CHECK(fd >= 0);
    staged = strdup(duckhts_output_path("directory", &owner));
    CHECK(staged != NULL);
    CHECK(close(fd) == 0);
    CHECK(duckhts_output_finish("directory", &owner, 1) < 0);
    CHECK(stat("directory", &st) == 0 && S_ISDIR(st.st_mode));
    CHECK(access(staged, F_OK) < 0 && errno == ENOENT);
    free(staged);

    fd = duckhts_output_open("identity", 0, &owner);
    CHECK(fd >= 0);
    CHECK(close(fd) == 0);
    CHECK(rename("identity", "moved") == 0);
    FILE *replacement = fopen("identity", "wb");
    CHECK(replacement != NULL);
    CHECK(fputs("unrelated", replacement) >= 0);
    CHECK(fclose(replacement) == 0);
    duckhts_output_cleanup("identity", &owner);
    CHECK(expect_contents("identity", "unrelated") == 0);

    CHECK(remove("output") == 0);
    CHECK(remove("identity") == 0);
    CHECK(remove("moved") == 0);
    CHECK(rmdir("directory") == 0);
    return 0;
}
