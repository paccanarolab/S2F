from Utils import ColourClass, FancyApp


class Measures(FancyApp.FancyApp):
    def __init__(self, args):
        super(Measures, self).__init__()
        self.colour = ColourClass.bcolors.OKGREEN
        self.config_file = args.config_file

    def run(self):
        pass
