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
        context = ET.iterparse(io.BytesIO(xmlstring.encode("utf-8")), events=('end',), tag=tag)
        texts = self.fast_xml_iter(context, lambda elem: None)
        return texts

    def fast_xml_iter(self, context, func, *args, **kwargs):
        """
        create a json tag per each document with the following structure
        {
            "datecreated": date_created,
            "article": {
                "abstract": abstract_text,
                "journal": {
                    "title": title
                    "volume": volume,
                    "issue": issue
                },
                "articletitle": article_title,
                "authorlist": [{"lastname":lastname,
                                "forename": forename,
                                "initials":initials},...]
            },
            "country": country
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
                pubmedarticledict["article"]["abstract"] = abstract_text


            date_tag = elem.findall(".//DateCreated")

            year, month, day = date_tag[0].getchildren()

            date_created = datetime(year=int(year.text), month=int(month.text), day=int(day.text))
            pubmedarticledict["datecreated"] = date_created

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

            article_title_tag = elem.findall(".//ArticleTitle")
            if article_title_tag and article_title_tag[0].text:
                article_title = article_title_tag[0].text
                pubmedarticledict["article"]["articletitle"] = article_title


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
    a = "<families><family><name>Simpson</name><members><name>Hommer</name><name>Marge</name><name>Bart</name></members></family><family><name>Griffin</name><members><name>Peter</name><name>Brian</name><name>Meg</name></members></family></families>"
    xmlread = XMLRead("")
    xmlread.get_all_texts(a, "name")

    #a = io.BytesIO("sa".encode("utf-8"))