import lxml.etree as ET
import io
from datetime import datetime
import os
import gzip
import logging
import calendar

class XMLRead(object):
    def __init__(self, dirname):
        self.dirname = dirname
        self.preprocessing = None
        self.month_to_number = {v: k for k, v in enumerate(calendar.month_abbr)}
        self.logger = logging.getLogger(__name__)

    def __iter__(self):
        """
        iteration on a directory containing xml files
        get article in each file that has been processed
        and converted to json
        """
        for fn in os.listdir(self.dirname):
            if os.path.isdir(os.path.join(self.dirname, fn)):
                continue
            xmlstring = ""
            try:
                if fn.endswith(".gz"):
                    xmlstring = gzip.open(os.path.join(self.dirname, fn)).read()
                elif fn.endswith(".xml"):
                    xmlstring = open(os.path.join(self.dirname, fn)).read()
                else:
                    logging.exception("Unknown type file(s)")
            except:
                logging.exception("Unknown type file(s)")
                continue
            article_data = self.get_all_texts(xmlstring, tag="PubmedArticle")
            for article in article_data:
                yield article

    def get_all_texts(self, xmlstring, tag):
        """
        gets all the text in a specified tag
        :param xmlstring:
        :param tag:
        :return:
        """
        if type(xmlstring) is bytes:
            xmlstring = xmlstring.decode('utf-8')

        try:
            context = ET.iterparse(io.BytesIO(xmlstring.encode("utf-8")), events=('end',), tag=tag)
        except:
            logging.exception("something is wrong with the xmlstring")
            raise

        texts = self.fast_xml_iter(context, lambda elem: None)
        return texts

    def fast_xml_iter(self, context, func, *args, **kwargs):
        """
        create a json tag per each document with the following structure
        {
        "PMID": PMID,
        "datecreated": date_created,
        "datecompleted": date_completed
        "article": {
            "abstract": abstract_text,
            "journal": {
                "title": title
                "volume": volume,
                "issue": issue,
                "ISSN": ISSN,
                "pubdate": pubdate
            },
            "articletitle": article_title,
            "authorlist": [{"lastname":lastname,
                            "forename": forename,
                            "initials":initials},...]
        },
        "country": country,
        "medlinejournalinfo": {
            "country": country,
            "medlineta": medline_ta,
            "ISSNlinking": ISSN_linking,
        },
        }
        :param context: contains events and elements tags
        :param func:
        :param args:
        :param kwargs:
        :return:
        """
        all_texts = []
        for event, elem in context:
            abstract_tag = elem.findall(".//AbstractText")
            pubmedarticledict = {}
            pubmedarticledict["article"] = {}
            pubmedarticledict["article"]["journal"] = {}
            if abstract_tag and abstract_tag[0].text:
                abstract_text = abstract_tag[0].text
                pubmedarticledict["article"]["abstract"] = abstract_text.lower()

            PMID_tag = elem.findall(".//PMID")
            if PMID_tag and PMID_tag[0].text:
                PMID = PMID_tag[0].text
                pubmedarticledict["PMID"] = PMID


            date_created_tag = elem.findall(".//DateCreated")
            if date_created_tag:
                year, month, day = date_created_tag[0].getchildren()
                date_created = datetime(year=int(year.text), month=int(month.text), day=int(day.text))
                pubmedarticledict["datecreated"] = date_created


            date_completed_tag = elem.findall(".//DateCompleted")
            if date_completed_tag:
                year, month, day = date_completed_tag[0].getchildren()
                date_completed = datetime(year=int(year.text), month=int(month.text), day=int(day.text))
                pubmedarticledict["datecompleted"] = date_completed


            title_tag = elem.findall(".//Title")
            if title_tag and title_tag[0].text:
                title = title_tag[0].text
                pubmedarticledict["article"]["journal"]["title"] = title

            volume_tag = elem.findall(".//Volume")
            if volume_tag and volume_tag[0].text:
                volume = volume_tag[0].text
                pubmedarticledict["article"]["journal"]["volume"] = volume

            issue_tag = elem.findall(".//Issue")
            if issue_tag and issue_tag[0].text:
                issue = issue_tag[0].text
                pubmedarticledict["article"]["journal"]["issue"] = issue


            ISSN_tag = elem.findall(".//ISSN")
            if ISSN_tag and ISSN_tag[0].text:
                ISSN = ISSN_tag[0].text
                pubmedarticledict["article"]["journal"]["ISSN"] = ISSN

            pubdate_tag = elem.findall(".//PubDate")
            if pubdate_tag and len(pubdate_tag[0].getchildren())==2:
                year, month = pubdate_tag[0].getchildren()
                try:
                    pubdate = datetime(year=int(year.text), month=self.month_to_number[month.text], day=1)
                except:
                    pubdate = datetime(year=int(year.text), month=1, day=1)
                pubmedarticledict["article"]["journal"]["pubdate"] = pubdate



            article_title_tag = elem.findall(".//ArticleTitle")
            if article_title_tag and article_title_tag[0].text:
                article_title = article_title_tag[0].text
                pubmedarticledict["article"]["articletitle"] = article_title.lower()


            author_list_tag = elem.find(".//AuthorList")
            authors = []
            if author_list_tag is not None:
                for author_tag in author_list_tag.findall(".//Author"):
                    lastname_tag = author_tag.find("LastName")
                    foretname_tag = author_tag.find("ForeName")
                    initials_tag = author_tag.find("Initials")
                    if lastname_tag is not None and foretname_tag is not None and initials_tag is not None:
                        authors.append({"lastname": lastname_tag.text,
                                        "forename": foretname_tag.text,
                                        "initials": initials_tag.text})
                pubmedarticledict["article"]["authorlist"] = authors

            country_tag = elem.findall(".//Country")
            if country_tag:
                country = country_tag[0].text
                pubmedarticledict["country"] = country

            medline_journal_info_tag = elem.find(".//MedlineJournalInfo")
            if medline_journal_info_tag is not None:
                pubmedarticledict["medlinejournalinfo"] = {}

            for medline_journal_info_item_tag in medline_journal_info_tag.getchildren():

                if medline_journal_info_item_tag.tag == "Country":
                    pubmedarticledict["medlinejournalinfo"]["country"] = medline_journal_info_item_tag.text

                if medline_journal_info_item_tag.tag == "MedlineTA":
                    pubmedarticledict["medlinejournalinfo"]["medlineta"] = medline_journal_info_item_tag.text

                if medline_journal_info_item_tag.tag == "ISSNLinking":
                    pubmedarticledict["medlinejournalinfo"]["ISSNlinking"] = medline_journal_info_item_tag.text


            all_texts.append(pubmedarticledict)

            func(elem, *args, **kwargs)
            # It's safe to call clear() here because no descendants will be accessed
            elem.clear()
            # Also eliminate now-empty references from the root node to elem
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]
        del context
        return all_texts


if __name__ == "__main__":
    #xmlread = XMLRead(dirname = "/media/rohola/09183702400/home/Desktop/ftp")
    xmlread = XMLRead(dirname = "/home/rohola/tmp/a")
    for item in xmlread:
        print(item)

    #a = io.BytesIO("sa".encode("utf-8"))