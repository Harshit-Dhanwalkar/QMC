#include "export_schema.h"

#include <time.h>

int export_timestamp_utc(char *buf, size_t size) {
  if (!buf || size < EXPORT_TIMESTAMP_SIZE) {
    return -1;
  }

  const time_t now = time(NULL);
  if (now == (time_t)-1) {
    return -1;
  }

  struct tm parts;
  if (!gmtime_r(&now, &parts)) {
    return -1;
  }

  if (strftime(buf, size, "%Y-%m-%dT%H:%M:%SZ", &parts) == 0) {
    return -1;
  }

  return 0;
}
