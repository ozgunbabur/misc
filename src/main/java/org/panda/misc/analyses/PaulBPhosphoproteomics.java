package org.panda.misc.analyses;

import org.panda.resource.UniProtSequence;
import org.panda.utility.CollectionUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class PaulBPhosphoproteomics
{
	public static void main(String[] args)
	{
	}

	static String getModifications(String markedPeptide, int startPos)
	{
		markedPeptide = markedPeptide.substring(2, markedPeptide.length() - 2);
		List<String> sites = new ArrayList<>();
		String[] t = markedPeptide.split("\\*");
		for (int i = 0; i < t.length; i++)
		{
			char lett = t[i].charAt(t[i].length() - 1);
			int pos = startPos + t[i].length() - 1;

			sites.add(lett + "" + pos);

			startPos += t[i].length();
		}

		return CollectionUtil.merge(sites, "|");
	}

}
