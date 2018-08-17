from Bio import Entrez


def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='20000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results



if __name__ == "__main__":
    results = dict(search("breast AND cancer"))
    print(results)