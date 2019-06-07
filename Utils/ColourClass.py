def coloured_string(colour, string):
    return colour + string + bcolors.ENDC

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    FAIL = '\033[1;31;40m'
    #WARNING = '\033[5;33;40m'
    WARNING = '\033[1;31m'
    ENDC = '\033[0m'
    RESET = "\033[m"
    BOLD = "\033[1m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    BOLD_RED = "\033[1;31m"
    BOLD_GREEN = "\033[1;32m"
    BOLD_YELLOW = "\033[1;33m"
    BOLD_BLUE = "\033[1;34m"
    BOLD_MAGENTA = "\033[1;35m"
    BOLD_CYAN = "\033[1;36m"
    BG_RED = "\033[41m"
    BG_GREEN = "\033[42m"
    BG_YELLOW = "\033[43m"
    BG_BLUE = "\033[44m"
    BG_MAGENTA = "\033[45m"
    BG_CYAN = "\033[46m"
    
    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


if __name__ == '__main__':
    from itertools import product
    line = ''
#    for i in range(256):
#        line += '\033[' + str(i) + 'm' + str(i) + bcolors.ENDC + '\t'
#        if i%10 == 0:
#            line += '\n'
    line = ''
    for a, b, c in product(range(256),repeat=3):

        line += '\033[' + str(a) + ';' + str(b) + ';' + str(c) + 'm' + str(a) + ';' + str(b) + ';' + str(c) + bcolors.ENDC + '\t'
        if c%10 == 0:
            line += '\n'
    print(line)
