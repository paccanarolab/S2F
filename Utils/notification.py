import tweepy
import json

def get_api(cfg):
    auth = tweepy.OAuthHandler(cfg['consumer_key'], cfg['consumer_secret'])
    auth.set_access_token(cfg['access_token'], cfg['access_token_secret'])
    return tweepy.API(auth)

def load_configuration(config_file):
    data = None
    with open(config_file) as config:
        data = json.load(config)

    return data
    # example configuration file
    #{ 
    #    "consumer_key"        : "VALUE",
    #    "consumer_secret"     : "VALUE",
    #    "access_token"        : "VALUE",
    #    "access_token_secret" : "VALUE" 
    #}

if __name__ == "__main__":
    import argparse
    aparser = argparse.ArgumentParser()
    required_arguments = aparser.add_argument_group('required arguments')
    required_arguments.add_argument('-t','--to', help='recipient of the tweet', required=True)
    required_arguments.add_argument('-m','--message', help='message to send', required=True)
    required_arguments.add_argument('-c','--config', help='configuration file', required=True)
    
    args = aparser.parse_args()

    api = get_api(load_configuration(args.config))
    message = args.message.replace(r'\n','\n')
    status = api.send_direct_message(screen_name=args.to,text=message) 
    #status = api.send_direct_message(screen_name=args.to,text="asd\nasd") 
    # Yes, tweet is called 'status' rather confusing

