package org.panda.misc;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.miner.IDFetcher;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class HGNCChEBIIDFetcher implements IDFetcher
{
	public Set<String> fetchID(BioPAXElement ele)
	{
		Set<String> set = new HashSet<>();

		if (ele instanceof SmallMoleculeReference)
		{
			for (Xref xref : ((SmallMoleculeReference) ele).getXref())
			{
				if (xref instanceof UnificationXref && xref.getDb() != null && xref.getDb().equalsIgnoreCase("ChEBI") &&
					xref.getId() != null)
				{
					set.add(xref.getId());
				}
			}
		}
		else if (ele instanceof ProteinReference)
		{
			for (Xref xref : ((ProteinReference) ele).getXref())
			{
				if (xref instanceof RelationshipXref && xref.getDb() != null &&
					xref.getDb().equalsIgnoreCase("HGNC Symbol") && xref.getId() != null)
				{
					set.add(xref.getId());
				}
			}
		}

		return set;
	}
}
