import requests


def getReqest(endpoint, printConsole=False):

    from pprint import pprint

    BASE = 'http://127.0.0.1:5000/'

    response = requests.get(BASE + endpoint)

    if printConsole != False:
        pprint(response.json())

    return response


def putRequest(endpoint, printConsole=False):
    BASE = 'http://127.0.0.1:5000/'
    response = requests.put(BASE + endpoint, data={'hello':'world'})
    return response


if __name__ == '__main__':

    endpoints = ['introns']

    get = [getReqest(i, printConsole=True) for i in endpoints]

