/************************************************************
 * SAXMzxmlHandler.cpp
 *
 * Premiere version janvier 2005
 * Patrick Lacasse
 * placasse@mat.ulaval.ca
 *
 * 3/11/2005 (Brendan MacLean): Use eXpat SAX parser, and create SAXSpectraHandler
 *
 * See http://sashimi.sourceforge.net/software_glossolalia.html for
 * mzXML schema information.
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

#include "stdafx.h"
#include "saxmzxmlhandler.h"

SAXMzxmlHandler::SAXMzxmlHandler( vector<scan>& _vS)
	: SAXSpectraHandler(_vS)
{
	m_bInMsLevel2 = false;
	m_bInPrecursorMz =  false;
	m_bInPeaks = false;
	//added this (true by default)
	m_bNetworkData = false;
}

SAXMzxmlHandler::~SAXMzxmlHandler()
{
}

void SAXMzxmlHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
	if(isElement("scan", el))
	{
		if((m_cidLevel = atoi(getAttrValue("msLevel", attr))) == 2)
		{
			m_bInMsLevel2 = true;

			reset();	// Clean up for the next scan

			m_scanNum = atoi(getAttrValue("num", attr));
			m_tId = m_scanNum;
			while(m_sId.find(m_tId) != m_sId.end())	{
				m_tId++;
			}
			m_sId.insert(m_tId);
			m_peaksCount = atoi(getAttrValue("peaksCount", attr));
		}
	}
	else if(isElement("peaks", el))
	{
		m_bInPeaks = true;
		m_bLowPrecision = (strcmp("64", getAttrValue("precision", attr)) != 0);
	}
	else if(isElement("precursorMz", el))
	{
		m_bInPrecursorMz = true;
		m_precursorCharge = atoi(getAttrValue("precursorCharge", attr));
	}
}

void SAXMzxmlHandler::endElement(const XML_Char *el)
{
	if(isElement("peaks", el))
	{
		processData();
		m_bInPeaks = false;
	}
	else if(isElement("precursorMz", el))
	{
		processData();
		m_bInPrecursorMz = false;
	}
	else if(isElement("scan", el) && m_bInMsLevel2 == true)
	{
		pushSpectrum();
		m_bInMsLevel2 = false;
	}
}

void SAXMzxmlHandler::characters(const XML_Char *s, int len)
{
	if ((m_bInPeaks && m_cidLevel == 2) ||
		(m_bInPrecursorMz))
	{
		m_strData.append(s, len);
	}
}

void SAXMzxmlHandler::processData()
{
	if( m_bInPeaks && m_cidLevel == 2)
	{
		pushPeaks();
	}
	else if (m_bInPrecursorMz)
	{
		m_precursorMz = atof(m_strData.data());
	}

	m_strData.clear();
}
