import ast
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1] / "src" / "mcp_ic_tool" / "mcp_server.py"


def _is_mcp_tool(function: ast.AsyncFunctionDef) -> bool:
    for decorator in function.decorator_list:
        if not isinstance(decorator, ast.Call):
            continue
        target = decorator.func
        if (
            isinstance(target, ast.Attribute)
            and target.attr == "tool"
            and isinstance(target.value, ast.Name)
            and target.value.id == "mcp"
        ):
            return True
    return False


def _contains_union_annotation(node: ast.AST) -> bool:
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
        return True
    if isinstance(node, ast.Name) and node.id in {"Optional", "Union"}:
        return True
    if isinstance(node, ast.Attribute) and node.attr in {"Optional", "Union"}:
        return True
    return any(_contains_union_annotation(child) for child in ast.iter_child_nodes(node))


def _tool_functions() -> dict[str, ast.AsyncFunctionDef]:
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, ast.AsyncFunctionDef) and _is_mcp_tool(node)
    }


def test_mcp_tool_parameters_avoid_union_schema_annotations():
    offenders = []
    for function in _tool_functions().values():
        for arg in (*function.args.args, *function.args.kwonlyargs):
            if arg.arg == "ctx" or arg.annotation is None:
                continue
            if _contains_union_annotation(arg.annotation):
                offenders.append(f"{function.name}.{arg.arg}: {ast.unparse(arg.annotation)}")

    assert not offenders, (
        "MCP tool argument schemas must have a top-level JSON schema type. "
        "Avoid Optional/Union annotations on @mcp.tool parameters; use a typed "
        "argument with default None instead.\n"
        + "\n".join(offenders)
    )


def test_multi_value_mcp_inputs_use_explicit_array_parameters():
    tools = _tool_functions()

    dos_args = {
        arg.arg: ast.unparse(arg.annotation)
        for arg in tools["submit_dos_calculation"].args.args
        if arg.annotation is not None
    }
    band_args = {
        arg.arg: ast.unparse(arg.annotation)
        for arg in tools["submit_band_structure_calculation"].args.args
        if arg.annotation is not None
    }
    md_args = {
        arg.arg: ast.unparse(arg.annotation)
        for arg in tools["submit_md_calculation"].args.args
        if arg.annotation is not None
    }

    assert dos_args["input_url"] == "str"
    assert dos_args["input_urls"] == "List[str]"
    assert band_args["input_url"] == "str"
    assert band_args["input_urls"] == "List[str]"
    assert md_args["temperature"] == "float"
    assert md_args["temperatures"] == "List[float]"
