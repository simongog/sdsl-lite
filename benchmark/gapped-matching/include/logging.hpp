#pragma once

#define ELPP_STL_LOGGING
#define ELPP_NO_DEFAULT_LOG_FILE

#ifndef NDEBUG
#define ELPP_DEBUG_ASSERT_FAILURE
#define ELPP_STACKTRACE_ON_CRASH
#endif

#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

struct log {
    inline static void load_default_config(bool print_to_stdout)
    {
        el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
        el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
        el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
        el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
        el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);

        el::Configurations c;
        c.setGlobally(el::ConfigurationType::Enabled, "true");
        c.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s} %level: %msg");
        c.setGlobally(el::ConfigurationType::ToFile, "false");
        if (print_to_stdout)
            c.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
        else
            c.setGlobally(el::ConfigurationType::ToStandardOutput, "false");
        el::Loggers::reconfigureAllLoggers(c);
    }

    inline static void start_log(int argc, const char** argv, bool print_to_stdout = true)
    {
        START_EASYLOGGINGPP(argc, argv);
        load_default_config(print_to_stdout);
    }
};