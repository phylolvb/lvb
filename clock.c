#include "Lvb.h"
#include "clock.h"

void log_Time()
{
    time_t timer;
    char buffer[26];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(buffer, 26, "%H:%M (%d/%m/%Y)", tm_info);
    puts(buffer);

}

void logstim(void)
/* log start time with message */
{
    time_t tim;	/* time */

    tim = time(NULL);
    printf("Starting at: %s\n", ctime(&tim));

} /* end logstim() */
