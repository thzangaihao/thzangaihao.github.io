/* ============================================================
   文件名: global.js
   描述: 全站通用脚本
   模块:
   1. 配置 (MathJax)
   2. 工具 (复制功能)
   3. 逻辑 (底栏、Tab切换、哈希跳转)
   4. 插件 (Prism 代码高亮、代码块增强)
   5. 初始化 (DOMContentLoaded)
   ============================================================ */


/* ------------------------------------------------------------
   1. 全局配置 (Configuration)
   ------------------------------------------------------------ */
window.MathJax = {
    tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        displayMath: [['$$', '$$'], ['\\[', '\\]']]
    },
    options: {
        ignoreHtmlClass: 'tex2jax_ignore',
        processHtmlClass: 'tex2jax_process'
    }
};


/* ------------------------------------------------------------
   2. 全局工具函数 (Utilities)
   ------------------------------------------------------------ */

/* --- 复制状态管理变量 --- */
let copyFeedbackTimer = null;
let currentActiveElement = null;
let originalContentCache = "";

/**
 * 2.1 通用复制函数 (用于底栏邮箱等文本链接)
 * 特点: 点击后显示 "已复制" 文字反馈
 */
function copyToClipboard(text, element, event) {
    if (event) event.preventDefault(); // 阻止链接默认跳转

    // 如果点击了新元素，立即重置上一个元素的状态
    if (currentActiveElement && currentActiveElement !== element) {
        resetCopyFeedback();
    }
    // 如果连续点击，重置定时器
    if (copyFeedbackTimer) clearTimeout(copyFeedbackTimer);

    navigator.clipboard.writeText(text).then(() => {
        // 首次点击，缓存原始内容
        if (currentActiveElement !== element) {
            currentActiveElement = element;
            originalContentCache = element.innerHTML;
        }

        // 显示反馈
        element.innerHTML = '<i class="fas fa-check" style="color: #27ae60;"></i> 已复制';
        
        // 1.5秒后自动恢复
        copyFeedbackTimer = setTimeout(() => {
            resetCopyFeedback();
        }, 1500);

    }).catch(err => console.error('复制失败:', err));
}

/**
 * 内部辅助: 恢复元素原状
 */
function resetCopyFeedback() {
    if (currentActiveElement) {
        currentActiveElement.innerHTML = originalContentCache;
        currentActiveElement = null;
        originalContentCache = "";
    }
    if (copyFeedbackTimer) {
        clearTimeout(copyFeedbackTimer);
        copyFeedbackTimer = null;
    }
}

/**
 * 2.2 代码块专用复制函数
 * 特点: 仅改变图标状态，不显示大段文字
 */
function copyCode(text, btnElement) {
    navigator.clipboard.writeText(text).then(() => {
        const originalHtml = btnElement.innerHTML;
        btnElement.innerHTML = '<i class="fas fa-check"></i>'; // 仅显示对勾
        btnElement.classList.add('copied');
        
        setTimeout(() => {
            btnElement.innerHTML = originalHtml;
            btnElement.classList.remove('copied');
        }, 1500);
    }).catch(err => {
        console.error('复制失败:', err);
        btnElement.innerText = 'Error';
    });
}


/* ------------------------------------------------------------
   3. 页面逻辑函数 (Page Logic)
   ------------------------------------------------------------ */

/**
 * 3.1 加载公共底栏 (Footer)
 */
function loadFooter() {
    const path = window.footerPath || './footer.html';
    
    fetch(path)
        .then(res => {
            if (!res.ok) throw new Error(`Footer not found at: ${path}`);
            return res.text();
        })
        .then(data => {
            const el = document.getElementById('footer-placeholder');
            if (el) el.innerHTML = data;
        })
        .catch(err => console.error('底栏加载错误:', err));
}

/**
 * 3.2 索引页 Tab 切换
 */
function showCollection(collectionId, btnElement) {
    // 切换文章列表显隐
    document.querySelectorAll('.article-list').forEach(el => el.classList.remove('active'));
    const target = document.getElementById(collectionId);
    if (target) target.classList.add('active');

    // 切换按钮激活状态
    document.querySelectorAll('.collection-btn').forEach(el => el.classList.remove('active'));
    if (btnElement) btnElement.classList.add('active');
}

/**
 * 3.3 处理 URL 哈希跳转 (#data -> 自动打开对应 Tab)
 * 包含“强制回顶”逻辑，解决浏览器自动滚动导致的布局错位
 */
function handleHashNavigation() {
    const hash = window.location.hash.substring(1);
    if (!hash) return;

    // 查找对应按钮
    const targetBtn = document.querySelector(`button[onclick*="'${hash}'"]`);
    if (targetBtn) {
        showCollection(hash, targetBtn);
        
        // 核心修复：强制滚回顶部，对抗浏览器的锚点定位
        // 执行两次以确保覆盖浏览器的默认行为
        window.scrollTo(0, 0);
        setTimeout(() => window.scrollTo(0, 0), 10);
    }
}


/* ------------------------------------------------------------
   4. 代码高亮系统 (Code Highlighting)
   ------------------------------------------------------------ */

/**
 * 4.1 动态加载 Prism.js (按需加载，不浪费流量)
 */
/* ------------------------------------------------------------
   4. 代码高亮系统 (Code Highlighting)
   ------------------------------------------------------------ */

/**
 * 4.1 动态加载 Prism.js (增强版：修复路径问题 + 自动去缩进)
 */
function loadPrismHighlighter() {
    // 1. 检查是否有代码块
    const codeBlocks = document.querySelectorAll('pre code');
    if (codeBlocks.length === 0) return;

    // --- [新增] 预处理：去除 HTML 源码缩进带来的多余空格 ---
    codeBlocks.forEach(code => {
        // 获取原始文本
        let text = code.textContent;
        // 1. 去掉首尾空行
        text = text.replace(/^\n+|\n+$/g, '');
        // 2. 计算最小缩进量
        const lines = text.split('\n');
        const indent = lines.reduce((min, line) => {
            if (line.trim().length === 0) return min; // 跳过空行
            const spaces = line.match(/^ */)[0].length;
            return spaces < min ? spaces : min;
        }, Infinity);
        // 3. 每一行都切掉这个缩进
        if (indent > 0 && indent !== Infinity) {
            text = lines.map(line => line.substring(indent)).join('\n');
        }
        code.textContent = text;
    });
    // -------------------------------------------------------

    const cdnBase = 'https://cdnjs.cloudflare.com/ajax/libs/prism/1.29.0';
    const cssUrl = `${cdnBase}/themes/prism-tomorrow.min.css`;
    const coreUrl = `${cdnBase}/prism.min.js`;
    const loaderUrl = `${cdnBase}/plugins/autoloader/prism-autoloader.min.js`;

    // 2. 加载 CSS
    const link = document.createElement('link');
    link.rel = 'stylesheet';
    link.href = cssUrl;
    document.head.appendChild(link);

    // 3. 辅助函数：串行加载 JS
    const loadScript = (url) => new Promise((resolve, reject) => {
        const script = document.createElement('script');
        script.src = url;
        script.onload = resolve;
        script.onerror = reject;
        document.body.appendChild(script);
    });

    // 4. 按顺序加载：核心 -> 自动识别插件 -> [关键修复] 配置路径
    loadScript(coreUrl)
        .then(() => loadScript(loaderUrl))
        .then(() => {
            // [关键修复] 显式告诉 Autoloader 去哪里下载语言包
            // 解决 "Loading..." 失败或不显示颜色的核心
            if (window.Prism && window.Prism.plugins && window.Prism.plugins.autoloader) {
                window.Prism.plugins.autoloader.languages_path = `${cdnBase}/components/`;
            }
            
            // 强制触发一次全局高亮
            window.Prism.highlightAll();
        })
        .catch(err => console.error("Prism load failed:", err));
}

/**
 * 4.2 代码块增强 (添加 Mac 风格顶部栏)
 */
function enhanceCodeBlocks() {
    document.querySelectorAll('pre code').forEach(code => {
        const pre = code.parentElement;
        // 防止重复处理
        if (pre.parentElement.classList.contains('code-wrapper')) return;

        // 提取语言名称
        let langName = 'TEXT';
        const langClass = Array.from(code.classList).find(c => c.startsWith('language-'));
        if (langClass) langName = langClass.replace('language-', '').toUpperCase();

        // 1. 创建容器
        const wrapper = document.createElement('div');
        wrapper.className = 'code-wrapper';

        // 2. 创建顶部栏
        const header = document.createElement('div');
        header.className = 'code-header';

        // 左侧：语言标签
        const langLabel = document.createElement('span');
        langLabel.className = 'lang-label';
        langLabel.innerText = langName;

        // 右侧：复制按钮
        const copyBtn = document.createElement('button');
        copyBtn.className = 'copy-code-btn';
        copyBtn.innerHTML = '<i class="far fa-copy"></i>';
        copyBtn.onclick = () => copyCode(pre.innerText, copyBtn);

        header.appendChild(langLabel);
        header.appendChild(copyBtn);

        // 3. 组装 DOM
        pre.parentNode.insertBefore(wrapper, pre);
        wrapper.appendChild(header);
        wrapper.appendChild(pre);
    });
}

/* ------------------------------------------------------------
   表格响应式增强 (自动添加滚动容器)
   ------------------------------------------------------------ */
function makeTablesResponsive() {
    const tables = document.querySelectorAll('.content table');
    tables.forEach(table => {
        if (table.parentElement.classList.contains('table-wrapper')) return;
        const wrapper = document.createElement('div');
        wrapper.className = 'table-wrapper';
        table.parentNode.insertBefore(wrapper, table);
        wrapper.appendChild(table);
    });
}

/* ------------------------------------------------------------
   5. 初始化入口 (Initialization)
   ------------------------------------------------------------ */
window.addEventListener('DOMContentLoaded', () => {
    
    // [重要] 禁用浏览器自动恢复滚动位置 (配合 handleHashNavigation 使用)
    if ('scrollRestoration' in history) {
        history.scrollRestoration = 'manual';
    }

    // A. 通用加载逻辑
    loadFooter();
    enhanceCodeBlocks();    // 先构建代码块外壳
    loadPrismHighlighter(); // 后加载高亮逻辑
    makeTablesResponsive();

    // B. 索引页专用逻辑 (Category Page)
    if (document.body.classList.contains('category-page')) {
        handleHashNavigation();
    }
});