from Utils import ColourClass
import time


class FancyApp(object):

    def __init__(self):
        self.colour = ColourClass.bcolors.BOLD_BLUE
        self.__verbose__ = True

    @staticmethod
    def warning_text(string):
        return ColourClass.bcolors.WARNING + string + ColourClass.bcolors.ENDC

    @staticmethod
    def yell(colour, handle, *args, **kwargs):
        print(colour + '[' + time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + str(handle) + ']: ' +
              ColourClass.bcolors.ENDC + " ".join(map(str, args)), **kwargs)

    def warning(self, *args, **kwargs):
        self.tell(ColourClass.bcolors.WARNING + " ".join(map(str, args)) + ColourClass.bcolors.ENDC, **kwargs)

    def tell(self, *args, **kwargs):
        if self.__verbose__:
            print(self.colour + '[' + time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + str(self.__class__.__name__) + ']: ' +
                  ColourClass.bcolors.ENDC + " ".join(map(str, args)), **kwargs)


if __name__ == '__main__':
    fa = FancyApp()
    fa.tell('hello world!')
    fa.warning('test')
