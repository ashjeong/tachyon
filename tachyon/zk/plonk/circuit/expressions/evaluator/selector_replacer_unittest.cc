#include "tachyon/zk/plonk/circuit/expressions/evaluator/selector_replacer.h"

#include <memory>

#include "tachyon/zk/plonk/circuit/expressions/evaluator/test/evaluator_test.h"
#include "tachyon/zk/plonk/circuit/expressions/expression_factory.h"

namespace tachyon::zk {

using Expr = std::unique_ptr<Expression<GF7>>;

class SelectorReplacerTest : public EvaluatorTest {};

TEST_F(SelectorReplacerTest, Constant) {
  GF7 value = GF7::Random();
  std::unique_ptr<Expression<GF7>> expr =
      ExpressionFactory<GF7>::Constant(value);
  EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
}

TEST_F(SelectorReplacerTest, Selector) {
  Expr expr = ExpressionFactory<GF7>::Selector(Selector::Simple(1));
  std::vector<std::unique_ptr<Expression<GF7>>> replacements;
  replacements.push_back(ExpressionFactory<GF7>::Selector(Selector::Simple(2)));
  replacements.push_back(
      ExpressionFactory<GF7>::Selector(Selector::Complex(3)));

  EXPECT_DEATH(expr->ReplaceSelectors(replacements, true), "");
  EXPECT_EQ(*expr->ReplaceSelectors(replacements, false), *replacements[1]);

  expr = ExpressionFactory<GF7>::Selector(Selector::Complex(1));
  EXPECT_EQ(*expr->ReplaceSelectors(replacements, false), *replacements[1]);
}

TEST_F(SelectorReplacerTest, Fixed) {
  struct {
    int32_t rotation;
    size_t column_index;
  } tests[] = {
      {1, 0},
      {2, 1},
  };

  for (const auto& test : tests) {
    FixedQuery query(1, Rotation(test.rotation),
                     FixedColumn(test.column_index));
    Expr expr = ExpressionFactory<GF7>::Fixed(query);
    EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
  }
}

TEST_F(SelectorReplacerTest, Advice) {
  struct {
    int32_t rotation;
    size_t column_index;
  } tests[] = {
      {6, 2},
      {7, 3},
  };

  for (const auto& test : tests) {
    AdviceQuery query(1, Rotation(test.rotation),
                      AdviceColumn(test.column_index, Phase(0)));
    Expr expr = ExpressionFactory<GF7>::Advice(query);
    EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
  }
}

TEST_F(SelectorReplacerTest, Instance) {
  struct {
    int32_t rotation;
    size_t column_index;
  } tests[] = {
      {1, 1},
      {2, 2},
  };

  for (const auto& test : tests) {
    InstanceQuery query(1, Rotation(test.rotation),
                        InstanceColumn(test.column_index));
    Expr expr = ExpressionFactory<GF7>::Instance(query);
    EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
  }
}

TEST_F(SelectorReplacerTest, Challenges) {
  for (size_t i = 0; i < challenges_.size(); ++i) {
    Expr expr = ExpressionFactory<GF7>::Challenge(Challenge(i, Phase(0)));
    EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
  }
}

TEST_F(SelectorReplacerTest, Negated) {
  GF7 value = GF7::Random();
  Expr expr =
      ExpressionFactory<GF7>::Negated(ExpressionFactory<GF7>::Constant(value));
  EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
}

TEST_F(SelectorReplacerTest, Sum) {
  GF7 a = GF7::Random();
  GF7 b = GF7::Random();
  Expr expr = ExpressionFactory<GF7>::Sum(ExpressionFactory<GF7>::Constant(a),
                                          ExpressionFactory<GF7>::Constant(b));
  EXPECT_FALSE(expr->ContainsSimpleSelector());
}

TEST_F(SelectorReplacerTest, Product) {
  GF7 a = GF7::Random();
  GF7 b = GF7::Random();
  Expr expr = ExpressionFactory<GF7>::Product(
      ExpressionFactory<GF7>::Constant(a), ExpressionFactory<GF7>::Constant(b));
  EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
}

TEST_F(SelectorReplacerTest, Scaled) {
  GF7 a = GF7::Random();
  GF7 b = GF7::Random();
  Expr expr =
      ExpressionFactory<GF7>::Scaled(ExpressionFactory<GF7>::Constant(a), b);
  EXPECT_EQ(*expr, *expr->ReplaceSelectors({}, false));
}

}  // namespace tachyon::zk
