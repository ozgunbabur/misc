package org.panda.misc;

import org.biopax.paxtools.pattern.Constraint;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.miner.ConfigurableIDFetcher;
import org.biopax.paxtools.pattern.miner.ControlsStateChangeOfMiner;
import org.biopax.paxtools.pattern.miner.CustomFormat;
import org.biopax.paxtools.pattern.miner.SIFSearcher;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;

/**
 * @author Ozgun Babur
 */
public class SIFWithoutGenerics
{
	public static void main(String[] args) throws FileNotFoundException
	{
		String dir = "/home/babur/Documents/DARPA/BigMech/";
		SimpleIOHandler io = new SimpleIOHandler();
		Model model = io.convertFromOWL(new FileInputStream(dir + "REACH.owl"));

		SIFSearcher searcher = new SIFSearcher(
			new ConfigurableIDFetcher().seqDbStartsWithOrEquals("HGNC Symbol"),
			new MyCSCO());

		searcher.searchSIF(model, new FileOutputStream(dir + "REACH-wo-generics.sif"),
			new CustomFormat("mediator", "comments"));
	}

	public static Constraint equal(boolean equal)
	{
		return new Equality(equal);
	}


	static class MyCSCO extends ControlsStateChangeOfMiner
	{
		@Override
		public Pattern constructPattern()
		{
			Pattern p = super.constructPattern();
			p.add(equal(true), "controller ER", "generic controller ER");
			p.add(equal(true), "controller simple PE", "controller PE");
			p.add(equal(true), "input PE", "input simple PE");
			p.add(equal(true), "output PE", "output simple PE");
			p.add(equal(true), "changed generic ER", "changed ER");
			return p;
		}
	}
}
