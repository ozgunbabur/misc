package org.panda.misc;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.EntityReference;
import org.biopax.paxtools.model.level3.ModificationFeature;
import org.biopax.paxtools.model.level3.SequenceModificationVocabulary;
import org.biopax.paxtools.model.level3.Xref;
import org.panda.utility.TermCounter;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class InspectBioPAXFile
{
	public static void main(String[] args) throws FileNotFoundException
	{
		Set<String> genes = new HashSet<>();
		SimpleIOHandler io = new SimpleIOHandler();
		Model model = io.convertFromOWL(new FileInputStream("/home/babur/Documents/PC/PathwayCommons9.Detailed.BIOPAX.owl"));

		TermCounter tc = new TermCounter();
		for (ModificationFeature mf : model.getObjects(ModificationFeature.class))
		{
			SequenceModificationVocabulary type = mf.getModificationType();
			if (type != null)
			{
				for (String term : type.getTerm())
				{
					if (!term.startsWith("MOD_RES"))
					{
						tc.addTerm(term);
						if (term.contains("methyl") && term.contains("lysine"))
						{
							EntityReference er = mf.getEntityFeatureOf();
							if (er != null)
							{
								for (Xref xref : er.getXref())
								{
									if (xref.getDb() != null && xref.getDb().equalsIgnoreCase("hgnc symbol"))
									{
										genes.add(xref.getId());
									}
								}
							}
						}
					}
				}
			}
		}

		tc.print();
		System.out.println();
		genes.stream().sorted().forEach(System.out::println);
	}
}
