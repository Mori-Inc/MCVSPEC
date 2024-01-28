#include <XSUtil/Numerics/MathOperator.h>
#include <XSUtil/Utils/IosHolder.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSUtil/Utils/XSutility.h>
#include <cctype>
#include <cmath>
#include <stack>
#include <utility>

// MdefExpression
#include <XSFunctions/Utilities/MdefExpression.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/XSCall.h>
#include <XSFunctions/Utilities/XSModelFunction.h>

string MdefElementString[] = {"ENG", "ENGC", "NUM", "PARAM", "OPER", "UFUNC", "BFUNC", 
			      "LPAREN", "RPAREN", "COMMA", "XSMODEL", "CONXSMODEL",
			      "TABLEMODEL"};

// Access to the list of models

// Class MdefExpression::MdefExpressionError 

MdefExpression::MdefExpressionError::MdefExpressionError (const string& errMsg)
   :YellowAlert("\nMdefine Expression Error: ")
{
  *IosHolder::errHolder() << errMsg << std::endl;
}


// Class MdefExpression 
const string MdefExpression::s_allValidChars = string("_#,.+-/*^(){}:abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 \t\r\n");
MdefExpression::MathOpContainer MdefExpression::s_operatorsMap;
std::map<string,int> MdefExpression::s_precedenceMap;

MdefExpression::MdefExpression(const MdefExpression &right)
   : AbstractExpression(right),
     m_distinctParNames(right.m_distinctParNames),
     m_paramsToGet(right.m_paramsToGet),
     m_paramTokenIndex(right.m_paramTokenIndex),
     m_numericalConsts(right.m_numericalConsts),
     m_operators(right.m_operators),
     m_postfixElems(right.m_postfixElems),
     m_infixElems(right.m_infixElems),
     m_eLow(right.m_eLow),
     m_eHigh(right.m_eHigh),
     m_compType(right.m_compType),
     m_usingOtherMdefs(right.m_usingOtherMdefs),
     m_mdefName(right.m_mdefName),
     m_callsSpecDependentFunctions(right.m_callsSpecDependentFunctions)
{
   if (s_operatorsMap.empty())
      buildOperatorsMap();
}

MdefExpression::MdefExpression (std::pair<Real,Real> eLimits, const string& compType, const string& mdefName)
   : AbstractExpression(),
     m_distinctParNames(),
     m_paramsToGet(),
     m_paramTokenIndex(),
     m_numericalConsts(),
     m_operators(),
     m_postfixElems(),
     m_infixElems(),
     m_eLow(eLimits.first),
     m_eHigh(eLimits.second),
     m_compType(compType),
     m_usingOtherMdefs(),
     m_mdefName(mdefName),
     m_callsSpecDependentFunctions(false)
{
   if (s_operatorsMap.empty())
      buildOperatorsMap();
}


MdefExpression::~MdefExpression()
{
}


MdefExpression & MdefExpression::operator=(const MdefExpression &right)
{
   if (this != &right)
   {
      MdefExpression tmp(right);
      Swap(tmp);
   }
   return *this;
}


void MdefExpression::init (const string& exprString, bool removeWhitespace)
{
   AbstractExpression::init(exprString, removeWhitespace);
   convertForTableModels();
   convertToInfix();
   convertToPostfix();
}

void MdefExpression::Swap (MdefExpression& right)
{
   AbstractExpression::Swap(right);
   std::swap(m_distinctParNames,right.m_distinctParNames);
   std::swap(m_numericalConsts,right.m_numericalConsts);
   std::swap(m_postfixElems,right.m_postfixElems);
   std::swap(m_infixElems,right.m_infixElems);
   std::swap(m_operators,right.m_operators);
   std::swap(m_paramsToGet,right.m_paramsToGet);
   std::swap(m_paramTokenIndex,right.m_paramTokenIndex);
   std::swap(m_eLow,right.m_eLow);
   std::swap(m_eHigh,right.m_eHigh);
   std::swap(m_compType,right.m_compType);
   std::swap(m_usingOtherMdefs,right.m_usingOtherMdefs);
   std::swap(m_mdefName,right.m_mdefName);
   std::swap(m_callsSpecDependentFunctions,right.m_callsSpecDependentFunctions);
}

MdefExpression* MdefExpression::clone () const
{
   return new MdefExpression(*this);
}

const string& MdefExpression::allValidChars () const
{
   return s_allValidChars;
}

void MdefExpression::buildOperatorsMap ()
{
   using namespace Numerics;
   clearOperatorsMap();
   // s_operatorsMap will own the memory of these objects throughout
   // the program's lifetime.
   s_operatorsMap["+"] = new PlusOp();
   s_operatorsMap["-"] = new MinusOp();
   s_operatorsMap["*"] = new MultOp();
   s_operatorsMap["/"] = new DivideOp();
   s_operatorsMap["^"] = new PowOp();
   s_operatorsMap["max"] = new MaxOp();
   s_operatorsMap["min"] = new MinOp();
   s_operatorsMap["atan2"] = new Atan2Op();
   s_operatorsMap["@"] = new UnaryMinusOp();
   s_operatorsMap["exp"] = new ExpOp();
   s_operatorsMap["sin"] = new SinOp();
   s_operatorsMap["sind"] = new SinDOp();
   s_operatorsMap["cos"] = new CosOp();
   s_operatorsMap["cosd"] = new CosDOp();
   s_operatorsMap["tan"] = new TanOp();
   s_operatorsMap["tand"] = new TanDOp();
   s_operatorsMap["sinh"] = new SinhOp();
   s_operatorsMap["sinhd"] = new SinhDOp();
   s_operatorsMap["cosh"] = new CoshOp();
   s_operatorsMap["coshd"] = new CoshDOp();
   s_operatorsMap["tanh"] = new TanhOp();
   s_operatorsMap["tanhd"] = new TanhDOp();
   s_operatorsMap["log"] = new LogOp();
   s_operatorsMap["ln"] = new LnOp();
   s_operatorsMap["sqrt"] = new SqrtOp();
   s_operatorsMap["abs"] = new AbsOp();
   s_operatorsMap["int"] = new IntOp();
   s_operatorsMap["sign"] = new SignOp();
   s_operatorsMap["heaviside"] = new HOp();
   s_operatorsMap["boxcar"] = new BoxcarOp();
   s_operatorsMap["asin"] = new ASinOp();
   s_operatorsMap["acos"] = new ACosOp();
   s_operatorsMap["atan"] = new ATanOp();
   s_operatorsMap["asinh"] = new ASinhOp();
   s_operatorsMap["acosh"] = new ACoshOp();
   s_operatorsMap["atanh"] = new ATanhOp();
   s_operatorsMap["mean"] = new MeanOp();
   s_operatorsMap["dim"] = new DimOp();
   s_operatorsMap["smin"] = new SMinOp();
   s_operatorsMap["smax"] = new SMaxOp();
   s_operatorsMap["erf"] = new ErfOp();
   s_operatorsMap["erfc"] = new ErfcOp();
   s_operatorsMap["gamma"] = new GammaOp();
   s_operatorsMap["legendre2"] = new Legendre2Op();
   s_operatorsMap["legendre3"] = new Legendre3Op();
   s_operatorsMap["legendre4"] = new Legendre4Op();
   s_operatorsMap["legendre5"] = new Legendre5Op();

   s_precedenceMap["+"] = 0;
   s_precedenceMap["-"] = 0;
   s_precedenceMap["@"] = 0;
   s_precedenceMap["*"] = 1;
   s_precedenceMap["/"] = 1;
   s_precedenceMap["#"] = 1;
   s_precedenceMap["^"] = 2;
}

void MdefExpression::clearOperatorsMap ()
{
   MathOpContainer::iterator itOp = s_operatorsMap.begin();
   MathOpContainer::iterator itOpEnd = s_operatorsMap.end();
   while (itOp != itOpEnd)   
   {
      delete itOp->second;
      ++itOp;
   }
   s_operatorsMap.clear();
   s_precedenceMap.clear();
}

void MdefExpression::convertForTableModels ()
{
  // catch cases of table model definitions which will be of form atable{...}, mtable{...}
  // or etable{...}. at the beginning these will be in tokenList as multiple tokens
  // and this routine combines them into a single token

  const std::vector<AbstractExpression::TokenType>& oldTokenList = tokenList();
  std::vector<AbstractExpression::TokenType> newTokenList;
  size_t nTokens = oldTokenList.size();
  size_t iTokens(0);
  while ( iTokens < nTokens ) {
    AbstractExpression::TokenType curToken = oldTokenList[iTokens];
    if ( curToken.tokenString == "atable" || curToken.tokenString == "mtable" ||
	 curToken.tokenString == "etable" ) {
      AbstractExpression::Token tabletype = curToken.type;
      size_t tableloc = curToken.location;
      string tablestring = curToken.tokenString;
      iTokens++;
      curToken = oldTokenList[iTokens];
      if ( curToken.tokenString != "{" ) {
	// missing { throw error and exit
      }
      // now carry on through the tokens until we reach the }
      while ( curToken.tokenString != "}" && iTokens < nTokens ) {
	tablestring += curToken.tokenString;
	iTokens++;
	curToken = oldTokenList[iTokens];
      }
      if ( curToken.tokenString != "}" ) {
	// missing } so throw error and exit
      }
      tablestring += curToken.tokenString;
      AbstractExpression::TokenType tableToken(tabletype, tableloc, tablestring);
      newTokenList.push_back(tableToken);
    } else {
      newTokenList.push_back(curToken);
    }
    iTokens++;
  }

  tokenList(newTokenList);

}

void MdefExpression::convertToInfix ()
{
   // Basically convert tokens enumerated by AbstractExpression's Token
   // types into MdefExpression's ElementTypes enumerators, which are more
   // useful for analyzing equations.
   const std::vector<AbstractExpression::TokenType>& exprTokenList = tokenList();
   std::vector<size_t> nonNumberTokens;
   findTheNumbers(nonNumberTokens,m_numericalConsts);
   int prevIdx = -1;
   int idx = -2;
   size_t numCount = 0;
   m_paramTokenIndex.clear();
   for (size_t i=0; i<nonNumberTokens.size(); ++i)
   {
      idx = static_cast<int>(nonNumberTokens[i]);
      MdefExpression::ElementType mathType = ENG; // init is irrelevant.
      // Assume any sequence of skipped tokens have been bundled into
      // one and only one number.
      if (idx > prevIdx+1)
      {
         m_infixElems.push_back(NUM);
         mathType = NUM;
         ++numCount;
      }
      const TokenType& curTok = exprTokenList[idx];
      switch (curTok.type)
      {
         case WordExp:
            mathType = classifyWords(curTok.tokenString);
	    if ( mathType == PARAM ) {
	      m_paramTokenIndex.push_back(idx);
	    }
            break;
         case Lbrace:
            // Conditions for implied '*': '(' is preceded by
            // ')' OR Eng OR param name OR NUM.
            if (m_infixElems.size())
            {
	      MdefExpression::ElementType prev = m_infixElems[m_infixElems.size()-1];
               if (prev == RPAREN || prev == ENG || prev == PARAM 
                        || prev == NUM || prev == ENGC)
               {
                  m_infixElems.push_back(OPER);
                  m_operators.push_back("*");
               }
            }
            mathType = LPAREN;
            break;
         case Rbrace:
            // This may also include an implied '*', which is
            // inserted after exiting the switch block.
            mathType = RPAREN;
            break;
         case Plus:
            mathType = OPER;
            m_operators.push_back("+");
            break;
         case Minus:
            mathType = OPER;
            // Conditions for unary: First token OR preceded by
            // a *,/,(, or Comma.
            if (idx == 0)
               m_operators.push_back("@");
            else
            {
               Token preType = exprTokenList[idx-1].type;
               if (preType == Lbrace || preType == Star || 
                        preType == Slash || preType == Comma)
                  m_operators.push_back("@");
               else
                  m_operators.push_back("-");
            }
            break;
         case Star:
            mathType = OPER;
            m_operators.push_back("*");
            break;
         case Slash:
            mathType = OPER;
            m_operators.push_back("/");
            break;
         case Exp:
            mathType = OPER;
            m_operators.push_back("^");
            break;
         case Comma:
            mathType = COMMA;
            break;
         default:
            {
               string errMsg("Unrecognized symbol during infix parsing: ");
               errMsg += curTok.tokenString;
               throw MdefExpressionError(errMsg);
            }
            break;
      }
      m_infixElems.push_back(mathType);
      // Conditions for implied '*' after ')': ')' is followed by
      // a WordExp of any kind.  
      if (mathType == RPAREN && static_cast<int>(exprTokenList.size()) > idx+1)
      {
         if (exprTokenList[idx+1].type == WordExp)
         {
            m_infixElems.push_back(OPER);
            m_operators.push_back("*");
         }
      }
      prevIdx = idx;
   } // End non-NUM token loop

   // idx may be negative here, so we don't want to cast it into a size_t.
   if (idx+1 < static_cast<int>(exprTokenList.size()))
   {
      // Assume the last token must be a number.  (There should never
      // be more than one remaining number, a condition which should 
      // already have been filtered out by AbstractExpression parsing.)
      m_infixElems.push_back(NUM);
      ++numCount;
   }
   if (numCount != m_numericalConsts.size())
     throw MdefExpressionError("Last symbol is not a number");

   // now we need to do some rearrangement in the case of any CONXSMODEL instances. In general
   // a convolution model will be written before the things it acts on however we need to
   // change this so a convolution model will be written after the things it acts on.

   // Start 
   
   verifyInfix();

   std::ostringstream oss;
   oss << "Infix elements: ";
   for (size_t i=0; i<m_infixElems.size(); ++i)
     oss << MdefElementString[m_infixElems[i]] << " ";
   oss << std::endl << "Distinct parameter names: ";
   for (size_t i=0; i<m_distinctParNames.size(); ++i)
     oss << m_distinctParNames[i] << " ";
   oss << std::endl << "Parameters to get: ";
   for (size_t i=0; i<m_paramsToGet.size(); ++i)
     oss << m_paramsToGet[i] << " ";
   oss << std::endl << "Parameters token index: ";
   for (size_t i=0; i<m_paramTokenIndex.size(); ++i)
     oss << m_paramTokenIndex[i] << " ";
   oss << std::endl << "Numerical consts: ";
   for (size_t i=0; i<m_numericalConsts.size(); ++i)
     oss << m_numericalConsts[i] << " ";
   oss << std::endl << "Infix operators: ";
   for (size_t i=0; i<m_operators.size(); ++i)
     oss << m_operators[i] << " ";
   oss << std::endl;
   FunctionUtility::xsWrite(oss.str(), 40);
   
}



void MdefExpression::convertToPostfix ()
{
   // This will fill in the m_postfixElems vector and also rearrange 
   // m_operators into the order they will be called in postfix.
   // The vectors for pars and numbers need no such reordering since 
   // they are called in the same sequence for infix and postfix.
   using namespace std;
   using Numerics::MathOperator;
   // The int in the pair below is for storing operator precedence number.
   stack<pair<int, string> > opStack;
   vector<string> tmpOperators;
   const size_t nElems = m_infixElems.size();
   size_t opPos = 0;
   bool isPrevConXSModel = false;
   for (size_t i=0; i<nElems; ++i) {
     MdefExpression::ElementType curType = m_infixElems[i];
     switch (curType) {
     case OPER:
       {
	 string curOp(m_operators[opPos]);
         // If '*' operator involves an xspec convolution model,
         // change it to '#' symbol so that evaluate() function
         // will know it needs special treatment.
         if (curOp == string("*") && isPrevConXSModel)
         {
            curOp = "#";
            m_operators[opPos] = curOp;
            isPrevConXSModel = false;
         }
	 int prec = s_precedenceMap.find(curOp)->second;
	 if (!opStack.empty() && curOp != "^") {
	   pair<int,string> topOp = opStack.top();
	   int testPrec = topOp.first;
	   while (testPrec >= prec) {
	     m_postfixElems.push_back(OPER);
	     tmpOperators.push_back(topOp.second);
	     opStack.pop();
	     if (opStack.empty())
	       testPrec = -999; // cause immediate exit from loop
	     else {
	       topOp = opStack.top();
	       testPrec = topOp.first;
	     }
	   }               
	 }
	 opStack.push(make_pair(prec, curOp));
	 ++opPos;
       }
       break;
     case UFUNC:
     case BFUNC:
     case XSMODEL:
     case CONXSMODEL:
     case TABLEMODEL:
       // Unary, binary, and xspec function calls will be handled the same way.  
       // Treat them as an LPAREN that also happens to have an actual operator 
       // associated with it. This means giving it a precedence of -1 so that 
       // RPAREN thinks it's a matching LPAREN, BUT also give it the actual 
       // function name rather than the blank string.
       opStack.push(make_pair(-1,m_operators[opPos]));
       ++opPos;
       // Now skip the actual LPAREN that we know is following this.
       ++i;
       break;
     case LPAREN:
       // If in here, this is a stand-alone parenthesis not associated with a 
       // function call.  LPAREN is not an actual operator, but we still need 
       // to store it as a placeholder in the opStack.  Give it a dummy 
       // precedence that's lower than any real operator, and a blank string.
       opStack.push(make_pair(-1," "));
       break;
     case RPAREN:
       // Pop the operators stack until we come to the first left parenthesis,
       // which may or may not be tied to a function call.
       {
	 int testPrec = 0;
	 do {
	   pair<int,string> topOp = opStack.top();
	   testPrec = topOp.first;
	   if (topOp.second != " ") {
	     // This could be an LPAREN tied to a function call, if testPrec
	     // = -1.  From this point forward, BFUNC and UFUNC need only be 
	     // classified as an OPER.
	     m_postfixElems.push_back(OPER);
	     tmpOperators.push_back(topOp.second);
             if (testPrec == -1)
             {
                // The LPAREN is from a function call. Find out if
                // the function is an xspec convolution model, which
                // will need special treatment during the evaluation.
                try
                {
                   ComponentInfo compInfo=XSModelFunction::compMatchName(topOp.second);
                   if (compInfo.type()==string("con"))
                      isPrevConXSModel = true;
                }
                catch (YellowAlert&)
                {
                   // Will end up in here if operator is not an xspec function but
                   //  a math function.
                }
             }
	   }
	   opStack.pop();  
	 } while (testPrec != -1 && !opStack.empty());
       }
       break;
     case COMMA:
       // Similar to RPAREN case, but do NOT pop the function call LPAREN.
       {
	 pair<int,string> topOp = opStack.top();
	 int testPrec = topOp.first;
	 // We can safely assume the first LPAREN reached is the 
	 // crucial function call LPAREN.  Any intervening ones would
	 // already have been popped when their RPAREN was processed.
	 while (testPrec != -1) {
	   m_postfixElems.push_back(OPER);
	   tmpOperators.push_back(topOp.second);
	   opStack.pop();
	   topOp = opStack.top();
	   testPrec = topOp.first;
	 }
       }
       break;
     default:
       m_postfixElems.push_back(curType);
       break;
     }
   }

   // Any operators remaining on the stack must now be popped to output.
   while (!opStack.empty())
   {
      m_postfixElems.push_back(OPER);
      tmpOperators.push_back(opStack.top().second);
      opStack.pop();
   }

   // And now the costly copy...
   m_operators = tmpOperators;

   std::ostringstream oss;
   oss << "Postfix elements: ";
   for (size_t i=0; i<m_postfixElems.size(); ++i)
     oss << MdefElementString[m_postfixElems[i]] << " ";
   oss << std::endl << "Postfix operators: ";
   for (size_t i=0; i<m_operators.size(); ++i)
     oss << m_operators[i] << " ";
   oss << std::endl;
   FunctionUtility::xsWrite(oss.str(), 40);

}

MdefExpression::ElementType MdefExpression::classifyWords (const string& wordStr)
{
   // Can assume wordStr is not a number, those have been taken care of
   // already.  It must either be a function, 'e', 'E', (with optional '.'
   // if convolution) or parameter name (which must start with a letter). 

   MdefExpression::ElementType type = COMMA; // Init to something other than ENG/ENGC 
   const bool isCon = (m_compType == string("con")); 
   if (wordStr.size() == 1 && (wordStr[0] == 'e' || wordStr[0] == 'E'))
      type = isCon ? ENGC : ENG;       
   else if (isCon && (wordStr == string(".e") || wordStr == string(".E")))
      type = ENG; 

   if (type != ENG && type != ENGC)
   {  
      // Function match will be case-insensitive, but other than that it
      // must match exactly.
      string lcWord = XSutility::lowerCase(wordStr);
      MathOpContainer::const_iterator itFunc = s_operatorsMap.find(lcWord);
      if (itFunc !=  s_operatorsMap.end() && lcWord == itFunc->first)
      {
         type = (itFunc->second->nArgs() == 1) ? UFUNC : BFUNC;
         m_operators.push_back(lcWord);
      } else {
	// It's not a function so check for XSPEC model using full name resolution
	bool foundXSPECmodel(false);
	if ( wordStr.length() > 2 ) {
	  try {
	    if ( XSModelFunction::isExactMatchName(wordStr) ) {
	      ComponentInfo nameFound = XSModelFunction::compMatchName(wordStr);
	      if ( nameFound.isMdefineModel() ) {
		// store information for dependencies on other mdef components
		m_usingOtherMdefs.insert(XSutility::lowerCase(nameFound.name()));
	      }
	      if ( nameFound.isSpecDependent() ) {
		// this mdefine'd model must now be spectrum dependent
		m_callsSpecDependentFunctions = true;
	      }
	      if ( nameFound.type() == string("con") ) {
		type = CONXSMODEL;
	      } else {
		type = XSMODEL;
	      }
	      m_operators.push_back(nameFound.name());
	      foundXSPECmodel = true;
            } else {
	      // catch special case of table models (atable, mtable, or etable) where
	      // wordStr should include the filename specification as well
	      string testName = wordStr.substr(0,6);
	      if ( testName == "atable" || testName == "mtable" || testName == "etable" ) {
		type = TABLEMODEL;
		m_operators.push_back(wordStr);
		foundXSPECmodel = true;
	      }
	    }
          } catch(...) {
	  }
	}
	if ( !foundXSPECmodel ) {
	  // At this point, must assume wordStr is a parameter name.
	  if (!(isalpha(wordStr[0]) || wordStr[0] == '_' || 
		wordStr.find(":") != string::npos))
	    {
	      string errMsg("Illegal parameter name: ");
	      errMsg += wordStr;
	      throw MdefExpressionError(errMsg);
	    }
	  // Determine if this name has already been entered, and if so,
	  // what is its index.   This is O(N^2) (ugh), but can't 
	  // imagine N > 100.
	  size_t idx=0;
	  std::vector<string>::const_iterator itName = m_distinctParNames.begin();
	  std::vector<string>::const_iterator itNameEnd = m_distinctParNames.end();
	  while (itName != itNameEnd)
	    {
	      if (*itName == wordStr)
		break;
	      ++itName, ++idx;
	    }
	  if (itName == itNameEnd) {
	    m_distinctParNames.push_back(wordStr);
	  }
	  m_paramsToGet.push_back(idx);
	  type = PARAM;
	}
      } 
   } // end if not ENG/ENGC
   return type;
}

void MdefExpression::verifyInfix () const
{
   // AbstractExpression has already verified general-syntax Token
   // ordering, which is really most of the work.  This will perform 
   // some additional checking that is specific to equation requirements:
   // - All function and xspec model names must immediately be followed by a '('.
   // - Commas can only exist in binary function or xspec model calls
   // - The proper number of commas must exist in (the top level of)
   //   a function call (1 for binary, nPar-1 for xspec model).
   // - Comma cannot appear inside any pair of parentheses internal to
   //   the pair specifying the binary function call to which it belongs.

   const size_t nElems = m_infixElems.size();
   size_t funcCounter = 0; // Needed only for output messages.
   // The first rule is easy enough to check with a linear search.
   // While doing this, count up all the commas and function
   // calls.  That nCommas = nBinary + totalnPars - nXspecMods
   // is necessary but not sufficient for rules 2 and 3.  
   // But once rule 3 is verified independently, this becomes sufficient 
   // for proving rule 2.  Rule 4 is verified internally in verifyBFunc.
   size_t commaCount = 0;
   size_t binaryCount = 0;
   size_t totXsModCommas = 0;
   size_t nXspecMods = 0, nMdefMods = 0;
   std::vector<size_t> xsModCommas;
   for (size_t i=0; i<nElems; ++i)
   {
      MdefExpression::ElementType curType = m_infixElems[i];
      if (curType == UFUNC || curType == BFUNC)
      {
         if (i == nElems-1 || m_infixElems[i+1] != LPAREN)
         {
            string errMsg("A '(' must follow the call to: ");
            errMsg += m_operators[funcCounter];
            throw MdefExpressionError(errMsg);
         }
         ++funcCounter;
         if (curType == BFUNC)
            ++binaryCount;
      }
      else if (curType == XSMODEL || curType == CONXSMODEL || curType == TABLEMODEL)
      {
         if (i == nElems-1 || m_infixElems[i+1] != LPAREN)
         {
            string errMsg("A '(' must follow the call to: ");
	    errMsg += m_operators[funcCounter];
            throw MdefExpressionError(errMsg);
         }
         ++funcCounter;
         ++nXspecMods;
         string opName = m_operators[funcCounter-1];
	 size_t nModCommas;
	 if ( curType == TABLEMODEL ) {
	   int numberParams, numberSpectra, numberEnergies;
	   bool isAdditive, isRedshift, isEscale;
	   string filename = opName.substr(7,opName.length()-8);
	   int status = FunctionUtility::tableInfo(filename, numberParams, numberSpectra,
						   numberEnergies, isAdditive, isRedshift,
						   isEscale);
	   if ( status != 0 ) {
	     string errMsg = "Filename " + filename + " cannot be found.";
	     throw MdefExpressionError(errMsg);
	   }
	   if ( isRedshift ) numberParams++;
	   if ( isEscale ) numberParams++;
	   nModCommas = numberParams-1;
	 } else {
	   nModCommas = XSModelFunction::numberParameters(opName)-1;
	 }
         // Store the number of commas for each xspec model in left-to-right
         //   order.  We'll need this information in the verifyFuncCommas
         //   function.
         xsModCommas.push_back(nModCommas);
	 totXsModCommas += nModCommas;
      }
      else if (curType == OPER)
         ++funcCounter;
      else if (curType == COMMA)
         ++commaCount;
   }
   if (commaCount > binaryCount + totXsModCommas) {
     string errMsg = "Extra commas detected in expression.\n";
     errMsg += "A common cause is that an unknown model name is included in the expression\n";
     errMsg += "or a function has been given the wrong number of arguments.";
     throw AbstractExpressionError(errMsg);
   }

   if (binaryCount || nXspecMods || nMdefMods)
   {
      size_t idxElem = 0;
      size_t ixsModFunc = 0;
      while (idxElem < nElems)
      {
         // When verifyBFunc returns, idxElem will be set to the index
         // of the closing parenthesis to the function for which it was
         // originally called.
         if (m_infixElems[idxElem] == BFUNC || m_infixElems[idxElem] == XSMODEL
	     || m_infixElems[idxElem] == CONXSMODEL || m_infixElems[idxElem] == TABLEMODEL)
            verifyFuncCommas(&idxElem,&ixsModFunc,xsModCommas);
         ++idxElem;
      }
   }   
}

void MdefExpression::verifyFuncCommas (size_t* idxElem, size_t* ixsFunc,
          const std::vector<size_t>& xsModCommas) const
{
   // This is a recursive function that verifies that one and only
   // one comma is placed in a binary function call, or nPars-1
   // commas is placed in an Xspec model function call.  
      
   // It also verifies that the comma is not inside nested parentheses.
   // ASSUME idxElem initially is the index of a binary or model function,
   // and idxElem+1 is the index of an LPAREN.  These things should 
   // already have been verified.
   
   // When function returns, idxElem will be set to the closing
   // parenthesis of the binary or model function.
   size_t parenCount = 1;
   size_t commasFound = 0;
   size_t expectedCommas = 0;
   MdefExpression::ElementType funcType = m_infixElems[*idxElem];
   if (funcType == BFUNC)
   {
      expectedCommas = 1;
   }
   else if (funcType == XSMODEL || funcType == CONXSMODEL || funcType == TABLEMODEL)
   {
      expectedCommas = xsModCommas[*ixsFunc];
      ++(*ixsFunc);
   }
   else
   {
      throw RedAlert("Programming error: Wrong operator type in verifyFuncCommas()");
   }
   
   size_t idx = *idxElem+1;
   while (parenCount)
   {
      ++idx;
      // Can also assume that expression does not end with a '('.
      // Therefore it should always be safe to start 2 elements 
      // further along.
      if (idx >= m_infixElems.size())
         throw RedAlert("Programming error in MdefExpression::verifyFuncCommas()");

      MdefExpression::ElementType curType = m_infixElems[idx];
      if (curType == COMMA)
      {
         ++commasFound;
         if (commasFound > expectedCommas)
            throw MdefExpressionError("Function called with too many arguments");
         // This is the test which verifies rule 4 stated in verifyInfix.   
         if (parenCount != 1)
            throw MdefExpressionError("Misplaced comma in function call.");
      }
      else if (curType == BFUNC || curType == XSMODEL || curType == CONXSMODEL || curType == TABLEMODEL)
      {
         verifyFuncCommas(&idx,ixsFunc,xsModCommas);
      }
      else if (curType == LPAREN)
      {
         ++parenCount;
      }
      else if (curType == RPAREN)
      {
         --parenCount;
      }
   }
   if (commasFound < expectedCommas)
   {
      throw MdefExpressionError("Function called with too few arguments");
   }
   *idxElem = idx;
}


void MdefExpression::evaluate (const RealArray& energies, const RealArray& parameters, int spectrumNumber,
			       RealArray& flux, RealArray& fluxErr, const string& inInitString) const
{
  using Numerics::MathOperator;
  using namespace std;

  string initString = inInitString;
  if (energies.size() < 2) throw MdefExpressionError("Energy array must be at least size 2");
  if (m_compType == string("con")) {
     convolveEvaluate(energies, parameters, spectrumNumber, flux, fluxErr, initString);
     return;
  }

  const size_t nBins = energies.size() - 1;

  RealArray avgEngs(nBins);
  RealArray binWidths(nBins);
  for (size_t i=0; i<nBins; ++i) {
    avgEngs[i] = (energies[i+1]+energies[i])/2.0;
    binWidths[i] = fabs(energies[i+1]-energies[i]);
  }

  // A boolean flag is coupled to the resultsStack arrays to mark whether
  // or not the array includes a factor of 1/binWidth, arising from XS add
  // components.  The convolution operator needs to know about this.
  typedef pair<RealArray, bool> MarkedArray;
  stack<MarkedArray> resultsStack;

  stack<RealArray> xsConParVals;
  stack<XSCallBase*> xsConFunctions;
  int numPos = 0;
  int parPos = 0;
  int opPos = 0;

  // m_postfixElems will contain the following enum values:
  //    ENG, NUM, PARAM, OPER
  for (auto& curType : m_postfixElems) {

    switch (curType) {

    case ENG:
    // ENGC should never get in here, but if it does just treat
    // it like ENG.
    case ENGC:
      resultsStack.push(MarkedArray(avgEngs,false));
      break;

    case NUM:
      resultsStack.push(MarkedArray(RealArray(m_numericalConsts[numPos],nBins),false));
      numPos++;
      break;

    case PARAM:
      {
	Real parVal = parameters[m_paramsToGet[parPos]];
	resultsStack.push(MarkedArray(RealArray(parVal,nBins),false));
	parPos++;
      }
      break;

    case OPER:
      {
	string opName = m_operators[opPos];

	MathOpContainer::const_iterator itFunc;

	if ( (itFunc = s_operatorsMap.find(opName)) != s_operatorsMap.end() &&
	     opName == itFunc->first) {
	  // case that OPER is a math operator
	  const MathOperator& mathFunc = *(itFunc->second);
	  const size_t nArgs = mathFunc.nArgs();
	  if (nArgs == 1) {
	    if (resultsStack.empty()) {
	      string msg = "Trying to access empty stack in MdefExpression::evaluate()\n";
	      msg += "                Likely error in mdefine expression. Try using chatter 40 to check.";
	      throw YellowAlert(msg);
	    }
	    RealArray& top = resultsStack.top().first;
	    mathFunc(top);
	  } else if (nArgs == 2) {
	    // Note that second array is not a reference.
	    if (resultsStack.size() < 2) {
	      string msg = "Too few arguments in MdefExpression::evaluate()\n";
	      msg += "                Likely error in mdefine expression. Try using chatter 40 to check.";
	      throw YellowAlert(msg);
	    }
	    RealArray second = resultsStack.top().first;
	    const bool isSecondDividedByBinWidth = resultsStack.top().second;
	    resultsStack.pop();
	    RealArray& first = resultsStack.top().first;
	    const bool isFirstDividedByBinWidth = resultsStack.top().second;
	    resultsStack.top().second = (isFirstDividedByBinWidth || isSecondDividedByBinWidth);
	    mathFunc(first, second);
	  }
	} else if (opName == string("#")){
	  // Special case of xspec model function requiring a convolution operation.
	  if (xsConParVals.empty() || xsConFunctions.empty()) {
	    throw RedAlert("Programmer Error: Mdefine operation with Xspec convolution model has empty stack.");
	  }
	  RealArray& modFlux = resultsStack.top().first;
	  const bool isDividedByBinWidth = resultsStack.top().second;
	  RealArray modFluxErr;
	  RealArray params(xsConParVals.top());
	  xsConParVals.pop();
	  const XSCallBase& modFunc = *(xsConFunctions.top());
	  xsConFunctions.pop();
	  if (isDividedByBinWidth) modFlux *= binWidths;
	  modFunc(energies, params, spectrumNumber, modFlux, modFluxErr, initString);
	  if (isDividedByBinWidth) modFlux /= binWidths;
	} else if ( XSModelFunction::hasFunctionPointer(opName) ) {
	  // case that OPER is an xspec model
	  const XSCallBase& modFunc = *(XSModelFunction::functionPointer(opName));
	  // get the number of parameters
	  size_t nParams = XSModelFunction::numberParameters(opName);
	  RealArray params(nParams);
	  // parameters will be popped off the stack in reverse order
	  for (size_t iparam=0; iparam<nParams; iparam++) {
	    params[nParams-iparam-1] = resultsStack.top().first[0];
	    resultsStack.pop();
	  }
	  // find the type of this component
	  const ComponentInfo compInfo = XSModelFunction::compMatchName(opName);
	  const string sub_compType = compInfo.type();
	  RealArray modFlux, modFluxErr;
	  // If the component is an xspec conv model, do NOT call
	  // its function here.  Just store the par vals and
	  // function pointer for now.  The convolution will be performed
	  // in the '#' handler, when it has the necessary flux array
	  // to operate on.
	  if (sub_compType == string("con") && !compInfo.isMdefineModel()) {
	    xsConParVals.push(params);
	    xsConFunctions.push(XSModelFunction::functionPointer(opName));
	  } else {
	    bool dividedByBinWidths=false;
	    modFunc(energies, params, spectrumNumber, modFlux, modFluxErr, initString);
	    if ( !compInfo.isMdefineModel() ) {

	      // if not an mdef model and the component type is add, con or mix then divide by the bin width
	      if (sub_compType == string("add")) {
		modFlux /= binWidths;
		dividedByBinWidths = true;
	      }

	    } else {

	      // if an mdef model the component type is mul or pileup then divide by the bin width
	      // is this correct?
	      if (sub_compType != string("mul") && sub_compType != string("pileup") ) {
		modFlux /= binWidths;
		dividedByBinWidths = true;
	      }

	    }

	    // push the result on the stack
	    resultsStack.push(MarkedArray(modFlux,dividedByBinWidths));
	  }

	} else if ( opName.substr(0,6) == "atable" || opName.substr(0,6) == "mtable" ||
		    opName.substr(0,6) == "etable" ) {
	  // special case to handle opName being a table model
	  size_t lenString = opName.length();
	  string filename = opName.substr(7,lenString-8);
	  string tableType("add");
	  if ( opName.substr(0,1) == "m" ) tableType = "mul";
	  if ( opName.substr(0,1) == "e" ) tableType = "exp";
	  // we need the number of parameters
	  int numberParams, numberSpectra, numberEnergies;
	  bool isAdditive, isRedshift, isEscale;
	  int status = FunctionUtility::tableInfo(filename, numberParams, numberSpectra,
						  numberEnergies, isAdditive, isRedshift,
						  isEscale);
	  if ( status != 0 ) {
	    string errMsg = "Filename " + filename + " cannot be found.";
	    throw MdefExpressionError(errMsg);
	  }
	  if ( isRedshift ) numberParams++;
	  if ( isEscale ) numberParams++;
	  // pop the parameters of the stack in reverse order
	  RealArray params(numberParams);
	  for (size_t iparam=0; iparam<numberParams; iparam++) {
	    params[numberParams-iparam-1] = resultsStack.top().first[0];
	    resultsStack.pop();
	  }
	  RealArray modFlux, modFluxErr;
	  FunctionUtility::tableInterpolate(energies, params, filename, spectrumNumber,
					    modFlux, modFluxErr, initString, tableType, false);
	  bool dividedByBinWidths(false);
	  if ( tableType == "add" ) {
	    modFlux /= binWidths;
	    dividedByBinWidths = true;
	  }
	  // push the result on the stack
	  resultsStack.push(MarkedArray(modFlux,dividedByBinWidths));

	} else {
	  // No function with that name found!  Probably because the user deleted it
	  YellowAlert("Attempt to call unknown model "+opName+
		      ". Did you delete an MDEFINEd model with that name?");
	  resultsStack.push(MarkedArray(RealArray(0.0,nBins),false));
	}
	opPos++;
      }
      break;

    default:
      string errMsg("Programmer error: unrecognized element type in MdefExpression::evaluate. ");
      break;

    } // end of switch over token type
  } // end m_postfixElems loop

  if (resultsStack.size() != 1)
    throw RedAlert("Programmer error: MdefExpression::evaluate() stack should be of size 1 at end.");
  if (flux.size() != nBins)
    flux.resize(nBins);
  flux = resultsStack.top().first;

  if (m_compType == string("add")) {
    // Integrate over bin, assume val is constant across bin.
    flux *= binWidths;
  }
}

void MdefExpression::convolveEvaluate (const RealArray& energies, const RealArray& parameters,
				       int spectrumNumber, RealArray& flux, RealArray& fluxErr,
				       const string& initString) const
{
   using Numerics::MathOperator;

   const size_t nBins = energies.size() - 1;
   if (flux.size() != nBins)
      throw RedAlert("Flux array size mismatch in mdef convolve function.");

   // test whether we can use the singleConvolveEvaluate. this is true if the final operator
   // is an xspec convolution model and there are no energies or other xspec models in the
   // expression
   if ( isSingleConvolve() ) {
     singleConvolveEvaluate(energies, parameters, spectrumNumber, flux, fluxErr, initString);
     return;
   }

   RealArray avgEngs(nBins);
   RealArray binWidths(nBins);
   for (size_t i=0; i<nBins; ++i)
   {
      avgEngs[i] = (energies[i+1]+energies[i])/2.0;
      binWidths[i] = fabs(energies[i+1]-energies[i]);
   }

   // To avoid building the numerical constant and parameter arrays
   // nBin times, do it once and store outside the loop (unlike the
   // standard evaluate function).
   const size_t nConsts = m_numericalConsts.size();
   const size_t nParsToGet = m_paramsToGet.size();
   std::vector<RealArray> constArrays(nConsts);
   std::vector<RealArray> paramArrays(nParsToGet);
   for (size_t i=0; i<nConsts; ++i)
   {
      constArrays[i].resize(nBins, m_numericalConsts[i]);
   }
   for (size_t i=0; i<nParsToGet; ++i)
   {
      paramArrays[i].resize(nBins, parameters[m_paramsToGet[i]]);
   }

   RealArray convFlux(0.0,nBins);
   for (size_t iBin=0; iBin<nBins; ++iBin)
   {
     // If input and convolved fluxes are row vectors [....], then
     // convEngs is the matrix convEngs_ji = avgEngs_i - avgEngs_j.
     RealArray convEngs(avgEngs[iBin] - avgEngs);
     std::stack<RealArray> resultsStack;
     size_t numPos = 0;
     size_t parPos = 0;
     size_t opPos = 0;

     // m_postfixElems will contain the following enum values
     //    ENG, ENGC, NUM, PARAM, OPER
     for (size_t iElem=0; iElem<m_postfixElems.size(); iElem++) {

       MdefExpression::ElementType curType = m_postfixElems[iElem];

       switch(curType) {

       case ENG:
           resultsStack.push(avgEngs);
	   break;

       case ENGC:
           resultsStack.push(convEngs);
	   break;

       case NUM:
	 resultsStack.push(constArrays[numPos]);
	 ++numPos;
	 break;

       case PARAM:
	 resultsStack.push(paramArrays[parPos]);
	 ++parPos;
	 break;

       case OPER:
	 {
	   string opName = m_operators[opPos];
	   MathOpContainer::const_iterator itFunc;

	   if ( (itFunc = s_operatorsMap.find(opName)) != s_operatorsMap.end() && 
		opName == itFunc->first) {
	     // case that OPER is a math operator
	     const MathOperator& mathFunc = *(itFunc->second);
	     const size_t nArgs = mathFunc.nArgs();
	     if (nArgs == 1) {
	       if (resultsStack.empty()) {
                 throw RedAlert("Trying to access empty stack in MdefExpression::convolveEvaluate()");
	       }
	       RealArray& top = resultsStack.top();
	       mathFunc(top);
	     } else if (nArgs == 2) {
	       // Note that second array is not a reference.
	       if (resultsStack.size() < 2)
                 throw RedAlert("Programmer error: Too few args in MdefExpression::convolveEvaluate() stack");
	       RealArray second = resultsStack.top();
	       resultsStack.pop();
	       RealArray& first = resultsStack.top();
	       mathFunc(first, second);
	     }
	   } else if ( XSModelFunction::hasFunctionPointer(opName) ) {
	     // case that OPER is an xspec model
	     const XSCallBase& modFunc = *(XSModelFunction::functionPointer(opName));
	     // get the number of parameters
	     size_t nParams = XSModelFunction::numberParameters(opName);
	     RealArray params(nParams);
	     // parameters will be popped off the stack in reverse order
	     for (size_t iparam=0; iparam<nParams; iparam++) {
	       params[nParams-iparam-1] = resultsStack.top()[0];
	       resultsStack.pop();
	     }
	     const ComponentInfo compInfo = XSModelFunction::compMatchName(opName);
	     const string sub_compType = compInfo.type();
	     RealArray modFlux, modFluxErr;
	     if ( sub_compType == string("con") ) {
	       if ( resultsStack.size() > 0 ) {
		 modFlux = resultsStack.top();
		 modFlux *= binWidths;
		 resultsStack.pop();
	       } else {
		 throw YellowAlert("Attempt to use a convolution component with nothing to operate on.");
	       }
	     }

	     modFunc(energies, params, spectrumNumber, modFlux, modFluxErr, initString);
	     if ( !compInfo.isMdefineModel() ) {
	     
	       // if not an mdef model and the component type is add, con or mix then divide by the bin width
	       if ( sub_compType == string("con") || sub_compType == string("add") ||
		    sub_compType == string("mix") ) {
		 modFlux /= binWidths;
	       }

	     } else {

	       // if an mdef model the component type is mul or pileup then divide by the bin width
	       // is this correct?
	       if (sub_compType != string("mul") && sub_compType != string("pileup") ) {
		 modFlux /= binWidths;
	       }

	     }
	     // if the component type is add, con or mix then push the resulting
	     // flux divided by the bin width onto the stack else just push the
	     // result flux
	     if ( sub_compType == string("con") || sub_compType == string("add") ||
		  sub_compType == string("mix") ) {
	       modFlux /= binWidths;
	     }
	     resultsStack.push(modFlux);

	     // if this is a convolution model then want to jump the next operator
	     // and next postfixElem because the "*" is not necessary
	     if ( m_compType == string("con") ) {
	       ++opPos;
	       iElem++;
	     }
	   }
	   else // No function with that name found!  Probably because the user deleted it
	   {
	     YellowAlert("Attempt to call unknown model "+opName+
			 ". Did you delete an MDEFINEd model with that name?");
	     resultsStack.push(RealArray(0.0,nBins));
	   }
	   ++opPos;
	 }
	 break;

       default:
	 throw RedAlert("Programmer error: unrecognized element type in MdefExpression::convolveEvaluate.");
	 break;

       } // end of switch over token type

     } // end m_postfixElems loop

     if (resultsStack.size() != 1)
        throw RedAlert("Programmer error: MdefExpression::convolveEvaluate() stack should be of size 1 at end.");

     // fact is a column vector
     RealArray& fact = resultsStack.top();
     fact *= binWidths[iBin];
     // Now multiply row and col vectors for new flux.
     convFlux[iBin] = (flux*fact).sum();        

   } // end iBins loop

   flux = convFlux;
}

void MdefExpression::singleConvolveEvaluate (const RealArray& energies, const RealArray& parameters,
					     int spectrumNumber, RealArray& flux, RealArray& fluxErr,
					     const string& initString) const
{
   using Numerics::MathOperator;

   // this is for special case of an mdefine convolution model which consists only of a single
   // xspec convolution model (either built-in of mdefine'd) - useful for redefining model
   // parameters of a convolution model

   const size_t nBins = energies.size() - 1;
   if (flux.size() != nBins)
      throw RedAlert("Flux array size mismatch in mdef convolve function.");

   // The last entry in m_postfixElems must be the xspec convolution model
   // so we start by performing all the operations before that. Since these
   // do not involve energy arrays we just need a stack of RealArrays of size 1

   std::stack<RealArray> resultsStack;

   int numPos(0);
   int parPos(0);
   int opPos(0);

   for (size_t iElem=0; iElem<m_postfixElems.size()-1; ++iElem) {

     MdefExpression::ElementType curType = m_postfixElems[iElem];
     switch (curType) {

     case ENG:
     case ENGC:
       // these should not occur
       throw RedAlert("Programmer error in MdefExpression::singleConvolveEvaluate(): found ENG or ENGC.");
       break;

     case NUM:
       resultsStack.push(RealArray(m_numericalConsts[numPos],1));
       numPos++;
       break;

     case PARAM:
       {
	 Real parVal = parameters[m_paramsToGet[parPos]];
	 resultsStack.push(RealArray(parVal,1));
	 parPos++;
       }
       break;

     case OPER:
       {
	 string opName = m_operators[opPos];
	 // this better be in the math operator list
	 MathOpContainer::const_iterator itFunc;
	 if ( (itFunc = s_operatorsMap.find(opName)) != s_operatorsMap.end() &&
	      opName == itFunc->first ) {
	   const MathOperator& mathFunc = *(itFunc->second);
	   const size_t nArgs = mathFunc.nArgs();
	   if (nArgs == 1) {
	     if (resultsStack.empty()) {
	       string msg = "Trying to access empty stack in MdefExpression::singleConvolutionEvaluate()\n";
	       msg += "                Likely error in mdefine expression. Try using chatter 40 to check.";
	       throw YellowAlert(msg);
	     }
	     RealArray& top = resultsStack.top();
	     mathFunc(top);
	   } else if (nArgs == 2) {
	     if (resultsStack.size() < 2) {
	       string msg = "Too few arguments in MdefExpression::singleConvolutionEvaluate()\n";
	       msg += "                Likely error in mdefine expression. Try using chatter 40 to check.";
	       throw YellowAlert(msg);
	     }
	     // Note that second is not a reference because mathFunc returns its result in first
	     RealArray second = resultsStack.top();
	     resultsStack.pop();
	     RealArray& first = resultsStack.top();
	     mathFunc(first,second);
	   }
	 } else {
	   // this should not happen
	   string msg = "Programmer error in MdefExpression::singleConvolveEvaluate():";
	   msg +=  " OPER " + opName + " is not math function.";
	   throw RedAlert(msg);
	 }
	 opPos++;
       }
       break;

     default:
       // this should not happen
       string msg = "Programmer error in MdefExpression::singleConvolveEvaluate(): unknown element type ";
       msg += curType;
       throw RedAlert(msg);
       break;

     } // end of switch

   } // end of loop over elements

   // now do the final element which is the convolution operation
   // opPos is set up to point to it

   string opName = m_operators[opPos];

   // if this is not an xspec model something has gone wrong
   if ( !XSModelFunction::hasFunctionPointer(opName) ) {
     // this should not happen
     string msg = "Programmer error in MdefExpression::singleConvolveEvaluate().";
     msg += " OPER " + opName + " does not have an xspec function pointer.";
     throw RedAlert(msg);
   }

   // it is an xspec model
   const XSCallBase& modFunc = *(XSModelFunction::functionPointer(opName));
   // get the number of parameters
   size_t nParams = XSModelFunction::numberParameters(opName);
   RealArray params(nParams);
   // parameters will be popped off the stack in reverse order
   for (size_t iparam=0; iparam<nParams; iparam++) {
     params[nParams-iparam-1] = resultsStack.top()[0];
     resultsStack.pop();
   }

   modFunc(energies, params, spectrumNumber, flux, fluxErr, initString);

}

bool MdefExpression::isSingleConvolve() const
{
  // test whether the expression uses a single xspec model which is convolution.
  // it must have no ENG or ENGC operators. first loop over all the elements
  // except the last one

  size_t opPos = 0;
  for (size_t iElem=0; iElem<m_postfixElems.size()-1; iElem++) {

    MdefExpression::ElementType curType = m_postfixElems[iElem];

    if ( curType == ENG || curType == ENGC ) return false;

    if ( curType == OPER ) {
      string opName = m_operators[opPos];
      // if this is not a math operator then return false
      if ( s_operatorsMap.find(opName) == s_operatorsMap.end() ) return false;
    }

  }

  // now check that the final element is an xspec convolution model

  MdefExpression::ElementType curType = m_postfixElems[m_postfixElems.size()-1];
  if ( curType != OPER ) return false;

  string opName = m_operators[m_operators.size()-1];
  if ( !XSModelFunction::hasFunctionPointer(opName) ) return false;

  const ComponentInfo compInfo = XSModelFunction::compMatchName(opName);
  if ( compInfo.type() != string("con") ) return false;

  return true;

}

// Additional Declarations
