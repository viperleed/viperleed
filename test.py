import logging

class CustomLogFormatter(logging.Formatter):
    """Logging Formatter for level-dependent message formatting"""

    FORMATS = {
#        logging.DEBUG: "DBG: %(module)s: %(lineno)d: %(msg)s",
        logging.DEBUG: "dbg: %(msg)s",
        logging.INFO: "%(msg)s",
        logging.WARNING: "# WARNING: %(msg)s",
        logging.ERROR: "### ERROR ### in %(module)s:%(funcName)s:%(lineno)s\n"
                       "# %(msg)s \n#############",
        logging.CRITICAL: "### CRITICAL ### in %(module)s:%(funcName)s:"
                          "%(lineno)s\n# %(msg)s \n################",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
    

# #start logger, write to file:
# timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
# logname = 'test-'+timestamp+'.log'

# #create logger:
# logFormatter = CustomLogFormatter()
# logging.getLogger(__name__).setLevel(logging.INFO)
# #add console output to logger:
# consoleHandler = logging.StreamHandler()
# consoleHandler.setFormatter(logFormatter)
# logging.getLogger(__name__).addHandler(consoleHandler)
# #add file output to logger:
# fileHandler = logging.FileHandler(logname, mode="w")
# fileHandler.setFormatter(logFormatter)
# logging.getLogger(__name__).addHandler(fileHandler)
# logger = logging.getLogger(__name__)
# logger.info("Starting new log: "+logname+"\nTime of execution (UTC): "
#              +time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+"\n")


# rootlogger = logging.getLogger()
# rootlogger.setLevel(logging.DEBUG)
def main():
    logger = logging.getLogger(__name__)
    # logger.setLevel(logging.DEBUG)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger.setLevel(logging.DEBUG)
    
    logFormatter = CustomLogFormatter()
    
    fileHandler = logging.FileHandler('test.log')
    fileHandler.setFormatter(logFormatter)
    # fileHandler.setLevel(logging.INFO)
    
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    # consoleHandler.setLevel(logging.DEBUG)
    
    logger.addHandler(consoleHandler)
    logger.addHandler(fileHandler)
    
    # print(logger.handlers)
    logger.warning('information message 8')
    
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])
    
    # print(logger.handlers)
    logging.shutdown()

main()