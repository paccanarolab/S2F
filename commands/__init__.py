def predict(args):
    from commands.Predict import Predict
    p = Predict(args)
    p.run()


def install(args):
    from commands.Install import Install
    i = Install(args)
    i.run()


def measures(args):
    from commands.Measures import Measures
    m = Measures(args)
    m.run()


def combine(args):
    from commands.Combine import Combine
    c = Combine(args)
    c.run()


def diffuse(args):
    from commands.Diffuse import Diffuse
    d = Diffuse(args)
    d.run()


def seed_from_hmmer(args):
    from commands.SeedFromHMMER import SeedFromHMMER
    s = SeedFromHMMER(args)
    s.run()
